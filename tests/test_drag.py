"""F1: physics, history, integration and backward compatibility."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from kessler_sim import DragConfig, Simulation, BreakupParams, load_scenario
from kessler_sim import drag, orbits, protocol
from kessler_sim.state import ObjectTable

DAY = 86400.
YEAR = 365.25 * DAY


def simulation(altitudes=(400., 800.), am=.02, e=0., ndot=0., created=-np.inf):
    altitudes = np.asarray(altitudes)
    n = len(altitudes)
    z = np.zeros(n)
    el = orbits.Elements((altitudes + orbits.R_MEAN) / (1 - e), z + e, z + .7,
                         z + .2, z + .3, z + .4, z + ndot)
    table = ObjectTable()
    table.add(np.arange(1, n + 1), [f'object {j}' for j in range(n)], 'fragment',
              1., .1, el, 0., t_created=created, am_m2_per_kg=am)
    return Simulation(table)


def test_density_positive_continuous_and_decreasing():
    h = np.linspace(150, 2000, 4000)
    rho = drag.atmospheric_density(h)
    assert np.all(rho > 0) and np.all(np.diff(rho) < 0)
    assert drag.atmospheric_density(400) == pytest.approx(3.725e-12)
    assert drag.atmospheric_density(800) == pytest.approx(1.170e-14)
    assert drag.atmospheric_density(500-1e-6) == pytest.approx(drag.atmospheric_density(500+1e-6), rel=1e-6)


def test_area_to_mass_and_validation():
    np.testing.assert_allclose(drag.area_to_mass([.02, np.nan, 0], [1, 2, 1], [10, 100, 1]), [.02, np.pi*4/100, 0])
    with pytest.raises(ValueError):
        drag.area_to_mass(np.nan, 1, 0)
    with pytest.raises(ValueError):
        DragConfig(cd=-1)


def test_circular_analytic_rate_and_units():
    cfg = DragConfig()
    a = orbits.R_MEAN + 400
    am = .02
    da, de = drag.decay_rates(a, 0., am, cfg)
    # Independent circular result in SI, then convert m/s to km/s.
    expected = -cfg.cd * am * drag.atmospheric_density(400) * np.sqrt(orbits.MU*1e9*a*1000) / 1000
    assert da == pytest.approx(expected, rel=1e-12)
    assert de == 0
    assert drag.decay_rates(a, 0, am*2, cfg)[0] == pytest.approx(2*da)
    assert drag.decay_rates(a, 0, 0, cfg)[0] == 0


def test_eccentric_decay_matches_reference_ode():
    sim = simulation((350,), e=.08)
    cfg = DragConfig(max_axis_change_km=2.)
    initial = sim.table.elements.take([0])
    assert drag.decay_rates(initial.a[0], initial.e[0], .02, cfg)[1] < 0
    t = 90*DAY
    def rhs(_, y):
        return drag.decay_rates(y[0], y[1], .02, cfg)
    ref = solve_ivp(rhs, (0, t), [initial.a[0], initial.e[0]], rtol=1e-10, atol=1e-12)
    sim.enable_drag(cfg)
    sim.advance_drag(t)
    el = sim.table.elements_at(t)
    assert el.a[0] < initial.a[0] and 0 < el.e[0] < initial.e[0]
    assert el.a[0] == pytest.approx(ref.y[0, -1], abs=.15)
    assert el.e[0] == pytest.approx(ref.y[1, -1], abs=3e-5)


def test_lifetimes_and_protocol_reentry():
    sim = simulation()
    sim.enable_drag()
    sim.advance_drag(50*YEAR)
    t = sim.reentries[1]
    assert 30*DAY < t < 3*YEAR
    assert 2 not in sim.reentries
    assert sim.object_info(2, 50*YEAR)['orbit']['perigee_alt_km'] > 150
    assert sim.object_info(1, t-1)['alive']
    assert not sim.object_info(1, t)['alive']
    assert sim.object_info(1, t)['position_km'] is None
    assert sim.object_info(1, t)['orbit']['perigee_alt_km'] == pytest.approx(150, abs=1e-6)
    assert np.isnan(protocol.state_positions(sim, t)[0]).all()
    assert sim.take_changes().updated == [1]
    assert sim.table.has(1)  # history retained


def test_seek_and_segment_boundary_continuity():
    sim = simulation()
    t0 = 1000.
    old = sim.state(t0)
    sim.enable_drag(t=t0)
    first_end = sim.table.orbit_history.latest[0, drag.END]
    sim.advance_drag(200*DAY)
    np.testing.assert_allclose(sim.state(t0).r_km, old.r_km, atol=1e-8)
    r_before, _ = sim.table.states(first_end - .001, [0])
    r_after, _ = sim.table.states(first_end + .001, [0])
    assert np.linalg.norm(r_after-r_before) < .03
    assert sim.table.elements_at(t0).a[0] == 6771
    assert sim.state(t0-1).ids.tolist() == [1, 2]
    assert sim.drag_time == 200*DAY  # seek did not alter history


def test_rewind_reentry_and_recompute():
    sim = simulation()
    sim.enable_drag()
    sim.advance_drag(YEAR)
    death = sim.reentries[1]
    expected = sim.state(death-10).r_km.copy()
    sim.take_changes()
    sim.truncate_after(1000.)
    assert not sim.reentries and sim.table.is_alive(0, 1000)
    assert sim.take_changes().updated == [1]
    sim.advance_drag(YEAR)
    assert sim.reentries[1] == death
    np.testing.assert_array_equal(sim.state(death-10).r_km, expected)


def test_sampling_cadence_does_not_change_results():
    a, b = simulation(), simulation()
    a.enable_drag(); b.enable_drag()
    a.advance_drag(180*DAY)
    for day in range(1, 181):
        b.advance_drag(day*DAY)
    assert a.reentries == b.reentries
    for t in (0, 30*DAY, 100*DAY, 180*DAY):
        np.testing.assert_array_equal(a.state(t).r_km, b.state(t).r_km)
    assert a.table.orbit_history.segment_count == b.table.orbit_history.segment_count


def test_catalogue_ndot_replaced_and_outside_scope_unchanged():
    sim = simulation((400, 2500), am=0., ndot=1e-12)
    original = sim.state(1000)
    high_before = sim.state(30*DAY).r_km[1]
    sim.enable_drag(t=1000)
    sim.advance_drag(30*DAY)
    np.testing.assert_allclose(sim.state(1000).r_km, original.r_km, atol=1e-8)
    assert sim.table.elements_at(30*DAY).ndot[0] == 0
    assert sim.table.elements_at(30*DAY).a[0] == 6771
    np.testing.assert_allclose(sim.state(30*DAY).r_km[1], high_before, atol=1e-8)


def test_collision_fragments_and_history_rewrite():
    sim = load_scenario('iridium_cosmos_2009')
    event = sim.events[0]
    sim.enable_drag(t=event.t)
    sim.advance_drag(event.t+DAY)
    ids = sim.state(event.t+DAY).ids
    ev2 = sim.trigger_breakup(int(ids[0]), int(ids[1]), event.t+DAY,
                              BreakupParams(impact_speed_kms=10, l_min_m=.3))
    sim.advance_drag(event.t+2*DAY)
    for oid in ev2.created_ids:
        assert np.isfinite(sim.table.orbit_history.latest[sim.table.row(oid), drag.START]) or not sim.table.is_alive(sim.table.row(oid), event.t+2*DAY)
    sim.truncate_after(event.t+DAY-1)
    assert len(sim.events) == 1
    assert all(not sim.table.has(oid) for oid in ev2.created_ids)
    assert all(sim.table.is_alive(sim.table.row(oid), event.t+DAY-1) for oid in (ev2.id_a, ev2.id_b))
    sim.advance_drag(event.t+DAY)
    assert sim.table.orbit_history.latest.shape[0] == len(sim.table)
    # Rewrite before activation: discard drag configuration and all segments.
    sim.truncate_after(event.t-1)
    assert sim.drag_config is None and not sim.reentries
    assert sim.table.orbit_history.segment_count == 0
    assert len(sim.events) == 0


def test_future_creation_and_invalid_time():
    sim = simulation((400,), created=10*DAY)
    sim.enable_drag()
    sim.advance_drag(20*DAY)
    assert not sim.object_info(1, 9*DAY)['alive']
    assert sim.table.elements_at(10*DAY).a[0] == 6771
    assert sim.table.elements_at(20*DAY).a[0] < 6771
    with pytest.raises(ValueError):
        sim.advance_drag(0)
    with pytest.raises(ValueError):
        sim.trigger_breakup(1, 2, 21*DAY)
