import numpy as np
import pytest

from kessler_sim import orbits as ob


def _random_states(n, seed=0):
    rng = np.random.default_rng(seed)
    r = rng.normal(size=(n, 3))
    r *= (rng.uniform(6600, 42500, n) / np.linalg.norm(r, axis=1))[:, None]
    v = rng.normal(size=(n, 3))
    v_circ = np.sqrt(ob.MU / np.linalg.norm(r, axis=1))
    v *= (v_circ * rng.uniform(0.7, 1.3, n) / np.linalg.norm(v, axis=1))[:, None]
    return r, v


def _elements(a, e, i_deg, raan_deg=0.0, argp_deg=0.0, m_deg=0.0):
    return ob.Elements(np.atleast_1d(float(a)), np.atleast_1d(float(e)), np.radians(np.atleast_1d(i_deg)),
                       np.radians(np.atleast_1d(raan_deg)), np.radians(np.atleast_1d(argp_deg)),
                       np.radians(np.atleast_1d(m_deg)))


def test_state_elements_round_trip_random():
    r, v = _random_states(2000)
    el, bound = ob.state_to_elements(r, v)
    assert bound.all()
    r2, v2 = ob.elements_to_state(el, 0.0)
    assert np.abs(r2 - r).max() < 1e-6
    assert np.abs(v2 - v).max() < 1e-9


@pytest.mark.parametrize("r,v", [
    ((7000.0, 0.0, 0.0), (0.0, 7.546, 0.0)),            # circular, equatorial, prograde
    ((7000.0, 0.0, 0.0), (0.0, -7.546, 0.0)),           # circular, equatorial, retrograde
    ((7000.0, 0.0, 0.0), (0.0, 0.0, 7.546)),            # circular, polar
    ((0.0, 7000.0, 0.0), (-8.5, 0.0, 0.0)),             # elliptical, equatorial
    ((0.0, 0.0, 7000.0), (8.0, 1.0, 0.0)),              # over the pole
])
def test_round_trip_singular_geometries(r, v):
    el, bound = ob.state_to_elements(np.array([r]), np.array([v]))
    assert bound[0]
    r2, v2 = ob.elements_to_state(el, 0.0)
    assert np.allclose(r2[0], r, atol=1e-6)
    assert np.allclose(v2[0], v, atol=1e-9)


def test_unbound_state_is_flagged():
    _, bound = ob.state_to_elements(np.array([[7000.0, 0, 0]]), np.array([[0, 11.0, 0]]))
    assert not bound[0]


def test_two_body_orbit_closes_after_one_period():
    el = _elements(7500.0, 0.1, 51.6, 30.0, 40.0, 10.0)
    r0, v0 = ob.elements_to_state(el, 0.0, with_j2=False)
    r1, v1 = ob.elements_to_state(el, ob.period(el)[0], with_j2=False)
    assert np.allclose(r0, r1, atol=1e-6)
    assert np.allclose(v0, v1, atol=1e-9)


def test_energy_and_angular_momentum_constant_along_orbit():
    el = _elements(8000.0, 0.2, 63.0, 10.0, 20.0, 0.0)
    times = np.linspace(0, 86400, 50)
    r, v = ob.elements_to_state(ob.Elements.concat([el] * len(times)), times)
    energy = 0.5 * np.einsum("ij,ij->i", v, v) - ob.MU / np.linalg.norm(r, axis=1)
    h = np.linalg.norm(np.cross(r, v), axis=1)
    assert np.allclose(energy, -ob.MU / (2 * 8000.0), rtol=1e-10)
    assert np.allclose(h, h[0], rtol=1e-10)


def test_kepler_solver_high_eccentricity():
    m = np.linspace(0, 2 * np.pi, 1000)
    for e in (0.0, 0.5, 0.9, 0.97):
        ecc = ob.solve_kepler(m, np.full_like(m, e))
        assert np.abs(ecc - e * np.sin(ecc) - np.mod(m, 2 * np.pi)).max() < 1e-10


def test_j2_sun_synchronous_precession():
    # 700 km sun-synchronous orbit: RAAN advances ~0.9856 deg/day (once per year)
    el = _elements(ob.R_EQ + 700.0, 0.0, 98.19)
    raan_dot, _, _ = ob.secular_rates(el)
    assert np.degrees(raan_dot[0]) * 86400 == pytest.approx(360.0 / 365.2422, rel=0.01)


def test_j2_critical_inclination_freezes_perigee():
    el = _elements(26600.0, 0.74, np.degrees(np.arccos(np.sqrt(1 / 5))))   # Molniya-like, 63.43 deg
    _, argp_dot, _ = ob.secular_rates(el)
    assert abs(argp_dot[0]) < 1e-15


def test_j2_iss_nodal_regression():
    # ISS (~420 km, 51.6 deg): RAAN regresses ~5 deg/day
    el = _elements(ob.R_EQ + 420.0, 0.0005, 51.64)
    raan_dot, _, _ = ob.secular_rates(el)
    assert np.degrees(raan_dot[0]) * 86400 == pytest.approx(-5.0, abs=0.1)


def test_j2_moves_orbit_plane_but_keeps_radius():
    el = _elements(7000.0, 0.0, 98.0)
    r_no, _ = ob.elements_to_state(el, 86400.0, with_j2=False)
    r_j2, _ = ob.elements_to_state(el, 86400.0)
    assert np.linalg.norm(r_j2) == pytest.approx(7000.0, rel=1e-9)
    assert np.linalg.norm(r_j2 - r_no) > 50.0


def test_mean_motion_decay_shifts_along_track():
    el = _elements(6800.0, 0.0, 51.6)
    decayed = ob.Elements(el.a, el.e, el.i, el.raan, el.argp, el.m, ndot=np.array([1e-12]))
    dt = 2 * 86400.0
    r0, _ = ob.elements_to_state(el, dt)
    r1, _ = ob.elements_to_state(decayed, dt)
    expected_angle = 0.5 * 1e-12 * dt ** 2           # rad
    angle = np.arccos(np.clip(np.dot(r0[0], r1[0]) / (np.linalg.norm(r0) * np.linalg.norm(r1)), -1, 1))
    assert angle == pytest.approx(expected_angle, rel=1e-6)


def test_vectorised_time_argument():
    el = ob.Elements.concat([_elements(7000.0, 0.0, 0.0)] * 2)
    r, _ = ob.elements_to_state(el, np.array([0.0, ob.period(el)[0] / 4]), with_j2=False)
    assert np.allclose(r[0], [7000.0, 0, 0], atol=1e-6)
    assert np.allclose(r[1], [0, 7000.0, 0], atol=1e-6)


def test_suborbital_orbit_has_no_j2_drift():
    # near-vertical climb at 865 km: perigee deep inside the Earth -> plain two-body motion
    el, bound = ob.state_to_elements(np.array([[6266.56, 0.0, 3618.0]]), np.array([[1.535, 0.069, 1.341]]))
    assert bound[0] and ob.perigee_radius(el)[0] < ob.R_EQ
    for dt in (-300.0, -100.0, 0.0):
        r_j2, v_j2 = ob.elements_to_state(el, dt)
        r_2b, v_2b = ob.elements_to_state(el, dt, with_j2=False)
        assert np.allclose(r_j2, r_2b) and np.allclose(v_j2, v_2b)


def test_orbits_above_earth_still_get_j2():
    el = _elements(7000.0, 0.01, 98.0)
    raan_dot, _, _ = ob.secular_rates(el)
    assert abs(raan_dot[0]) > 0
