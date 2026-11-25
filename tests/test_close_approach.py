import json
import math
import random
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC / "collision_detection"))
sys.path.insert(0, str(SRC / "collision_physics"))
import close_approach as ca  # noqa: E402
import fragmentation as fr  # noqa: E402


# ---------------------------------------------------------------------------
# TCA i odległość minimalna (ruch liniowy w obrębie kroku)
# ---------------------------------------------------------------------------

def test_tca_head_on():
    t, d, v, d0 = ca.compute_tca_and_miss_distance(
        (-10.0, 0.0, 7000.0), (1.0, 0.0, 0.0),
        (10.0, 0.0, 7000.0), (-1.0, 0.0, 0.0), dt_window=20.0)
    assert t == pytest.approx(10.0)
    assert d == pytest.approx(0.0, abs=1e-9)
    assert v == pytest.approx(2.0)
    assert d0 == pytest.approx(20.0)


def test_tca_offset_inside_step():
    t, d, _, _ = ca.compute_tca_and_miss_distance(
        (-60.0, 1.2, 7100.0), (6.0, 0.0, 0.0),
        (60.0, 0.0, 7100.0), (-6.0, 0.0, 0.0), dt_window=20.0)
    assert t == pytest.approx(10.0)
    assert d == pytest.approx(1.2)


def test_tca_perpendicular_crossing():
    # A leci wzdłuż x, B wzdłuż y; przecięcie torów w (0,0,0) w różnych chwilach
    # dr(t) = (-10 + t, 5 + t - 15, 0)  ->  min przy t = 10, d = 0  (gdy B startuje z y=-10)
    t, d, v, _ = ca.compute_tca_and_miss_distance(
        (-10.0, 0.0, 0.0), (1.0, 0.0, 0.0),
        (0.0, -10.0, 0.0), (0.0, 1.0, 0.0), dt_window=100.0)
    # względnie: dr0 = (-10, 10, 0), dv = (1, -1, 0) -> spotkanie w t = 10
    assert t == pytest.approx(10.0)
    assert d == pytest.approx(0.0, abs=1e-9)
    assert v == pytest.approx(math.sqrt(2.0))


def test_tca_diverging_clamped_to_window_start():
    t, d, _, d0 = ca.compute_tca_and_miss_distance(
        (-10.0, 0.0, 0.0), (-1.0, 0.0, 0.0),
        (10.0, 0.0, 0.0), (1.0, 0.0, 0.0), dt_window=5.0)
    assert t == 0.0
    assert d == pytest.approx(d0) == pytest.approx(20.0)


def test_tca_closest_after_window_clamped_to_window_end():
    t, d, _, _ = ca.compute_tca_and_miss_distance(
        (-100.0, 0.0, 0.0), (1.0, 0.0, 0.0),
        (100.0, 0.0, 0.0), (-1.0, 0.0, 0.0), dt_window=5.0)
    assert t == pytest.approx(5.0)
    assert d == pytest.approx(190.0)


def test_tca_zero_relative_velocity():
    t, d, v, d0 = ca.compute_tca_and_miss_distance(
        (0.0, 0.0, 0.0), (7.0, 1.0, 0.0),
        (3.0, 4.0, 0.0), (7.0, 1.0, 0.0), dt_window=5.0)
    assert (t, v) == (0.0, 0.0)
    assert d == d0 == pytest.approx(5.0)


def test_tca_symmetric_in_object_order():
    rng = random.Random(3)
    for _ in range(50):
        ra = tuple(rng.uniform(-50, 50) for _ in range(3))
        rb = tuple(rng.uniform(-50, 50) for _ in range(3))
        va = tuple(rng.uniform(-8, 8) for _ in range(3))
        vb = tuple(rng.uniform(-8, 8) for _ in range(3))
        assert ca.compute_tca_and_miss_distance(ra, va, rb, vb, 5.0) == pytest.approx(
            ca.compute_tca_and_miss_distance(rb, vb, ra, va, 5.0))


def test_tca_matches_brute_force_sampling():
    rng = random.Random(11)
    for _ in range(50):
        ra = tuple(rng.uniform(-50, 50) for _ in range(3))
        rb = tuple(rng.uniform(-50, 50) for _ in range(3))
        va = tuple(rng.uniform(-8, 8) for _ in range(3))
        vb = tuple(rng.uniform(-8, 8) for _ in range(3))
        _, d, _, _ = ca.compute_tca_and_miss_distance(ra, va, rb, vb, 5.0)
        sampled = min(
            ca.magnitude(ca.sub(ca.add_scaled(ra, va, k * 0.005), ca.add_scaled(rb, vb, k * 0.005)))
            for k in range(1001)
        )
        assert d <= sampled + 1e-9
        assert d == pytest.approx(sampled, abs=1e-3)


# ---------------------------------------------------------------------------
# Energetyka i klasyfikacja
# ---------------------------------------------------------------------------

def test_energetics_known_values():
    # 1 kg w 1000 kg przy 10 km/s: EMR = 0.5 * 1 * (1e4)^2 / 1e6 g = 50 J/g
    ke_mj, emr, cata = ca.compute_collision_energetics(1000.0, 1.0, 10.0)
    assert emr == pytest.approx(50.0)
    assert cata
    # energia w układzie środka masy (masa zredukowana)
    assert ke_mj == pytest.approx(0.5 * (1000.0 / 1001.0) * 1e8 / 1e6)


def test_energetics_order_independent():
    assert ca.compute_collision_energetics(10.0, 900.0, 3.0) == ca.compute_collision_energetics(900.0, 10.0, 3.0)


def test_energetics_catastrophic_threshold_inclusive():
    # EMR = 40 J/g dokładnie: 0.5 * 1 kg * v^2 / 1000 g = 40  ->  v = sqrt(80000) m/s
    v_kms = math.sqrt(80000.0) / 1000.0
    _, emr, cata = ca.compute_collision_energetics(1.0, 1.0, v_kms)
    assert emr == pytest.approx(40.0)
    assert cata


@pytest.mark.parametrize("m1,m2", [(0.0, 100.0), (100.0, 0.0), (-1.0, 5.0)])
def test_energetics_non_positive_mass(m1, m2):
    assert ca.compute_collision_energetics(m1, m2, 10.0) == (0.0, 0.0, False)


@pytest.mark.parametrize("m1,m2,v", [
    (500.0, 500.0, 7.0), (1500.0, 1.5, 10.0), (900.0, 560.0, 11.7), (1000.0, 0.01, 10.0), (500.0, 1.5, 0.6),
])
def test_energetics_consistent_with_fragmentation_model(m1, m2, v):
    _, emr_ca, cata_ca = ca.compute_collision_energetics(m1, m2, v)
    _, emr_fr, cata_fr = fr.collision_energetics(m1, m2, v)
    assert emr_ca == pytest.approx(emr_fr)
    assert cata_ca == cata_fr


@pytest.mark.parametrize("d,hbr,expected", [
    (0.001, 0.003, "PHYSICAL_COLLISION"),
    (0.003, 0.003, "PHYSICAL_COLLISION"),
    (0.5, 0.003, "CRITICAL_CONJUNCTION"),
    (ca.CRITICAL_CONJUNCTION_KM, 0.003, "CRITICAL_CONJUNCTION"),
    (ca.CRITICAL_CONJUNCTION_KM + 1e-6, 0.003, "CLOSE_APPROACH"),
    (9.0, 0.003, "CLOSE_APPROACH"),
])
def test_classify_event(d, hbr, expected):
    assert ca.classify_event(d, hbr) == expected


@pytest.mark.parametrize("a,b,expected", [
    (("STARLINK-1007", "active"), ("STARLINK-30001", "active"), True),
    (("ONEWEB-0012", "active"), ("ONEWEB-0600", "active"), True),
    (("STARLINK-1007", "active"), ("ONEWEB-0600", "active"), False),
    (("STARLINK-1007", "active"), ("STARLINK DEB", "debris"), False),
    (("ISS (ZARYA)", "active"), ("CSS (TIANHE)", "active"), False),
    (("starlink-1", "active"), ("Starlink-2", "active"), True),
])
def test_same_constellation(a, b, expected):
    obj_a = {"name": a[0], "category": a[1]}
    obj_b = {"name": b[0], "category": b[1]}
    assert ca.is_same_constellation_formation(obj_a, obj_b) is expected


# ---------------------------------------------------------------------------
# Siatka przestrzenna
# ---------------------------------------------------------------------------

def test_grid_finds_all_pairs_closer_than_cell_size():
    rng = random.Random(5)
    cell = 100.0
    objs = [{"r": tuple(rng.uniform(6500, 7100) for _ in range(3))} for _ in range(400)]
    pairs = ca.screen_candidate_pairs(objs, cell)
    found = {frozenset(p) for p in pairs}
    for i in range(len(objs)):
        for j in range(i + 1, len(objs)):
            if ca.magnitude(ca.sub(objs[i]["r"], objs[j]["r"])) < cell:
                assert frozenset((i, j)) in found


def test_grid_no_duplicate_or_self_pairs():
    rng = random.Random(6)
    objs = [{"r": tuple(rng.uniform(-300, 300) for _ in range(3))} for _ in range(300)]
    pairs = ca.screen_candidate_pairs(objs, 100.0)
    assert all(i != j for i, j in pairs)
    assert len({frozenset(p) for p in pairs}) == len(pairs)


def test_grid_cell_boundary_pair_found():
    objs = [{"r": (ca.CELL_SIZE_KM - 0.1, 0.0, 7000.0)}, {"r": (ca.CELL_SIZE_KM + 0.1, 0.0, 7000.0)}]
    assert len(ca.screen_candidate_pairs(objs, ca.CELL_SIZE_KM)) == 1


def test_grid_negative_coordinates():
    objs = [{"r": (-0.1, -0.1, -0.1)}, {"r": (0.1, 0.1, 0.1)}]
    assert len(ca.screen_candidate_pairs(objs, 100.0)) == 1


def test_grid_cell_covers_worst_case_motion_in_one_step():
    # Para zgłaszana na początku kroku musi obejmować każdy obiekt, który w ciągu kroku
    # może zbliżyć się na próg: max v_rel (~15.5 km/s w LEO) * krok + próg < rozmiar komórki
    assert 15.5 * ca.PROPAGATION_STEP_SECONDS + ca.DEFAULT_SCREENING_THRESHOLD_KM < ca.CELL_SIZE_KM


# ---------------------------------------------------------------------------
# Propagacja RK2
# ---------------------------------------------------------------------------

def _circular(r_km=7000.0):
    v = math.sqrt(ca.MU_EARTH / r_km)
    period = 2 * math.pi * math.sqrt(r_km ** 3 / ca.MU_EARTH)
    return {"r": (r_km, 0.0, 0.0), "v": (0.0, v, 0.0)}, period


def _energy(o):
    return 0.5 * ca.dot(o["v"], o["v"]) - ca.MU_EARTH / ca.magnitude(o["r"])


def test_rk2_circular_orbit_closes_after_one_period():
    obj, period = _circular()
    n = round(period / ca.PROPAGATION_STEP_SECONDS)
    for _ in range(n):
        ca.propagate_step_rk2([obj], period / n)
    assert ca.magnitude(ca.sub(obj["r"], (7000.0, 0.0, 0.0))) < 2.0      # km
    assert ca.magnitude(obj["r"]) == pytest.approx(7000.0, abs=0.01)


def test_rk2_conserves_energy_and_angular_momentum_on_ellipse():
    obj = {"r": (6800.0, 0.0, 0.0), "v": (0.0, 7.9, 1.2)}
    e0 = _energy(obj)

    def h(o):
        r, v = o["r"], o["v"]
        return (r[1] * v[2] - r[2] * v[1], r[2] * v[0] - r[0] * v[2], r[0] * v[1] - r[1] * v[0])

    h0 = h(obj)
    for _ in range(1200):  # 100 min
        ca.propagate_step_rk2([obj], 5.0)
    assert abs((_energy(obj) - e0) / e0) < 1e-5
    assert ca.magnitude(ca.sub(h(obj), h0)) / ca.magnitude(h0) < 1e-5


def test_rk2_skips_objects_below_minimum_radius():
    obj = {"r": (6400.0, 0.0, 0.0), "v": (0.0, 7.0, 0.0)}
    ca.propagate_step_rk2([obj], 5.0)
    assert obj == {"r": (6400.0, 0.0, 0.0), "v": (0.0, 7.0, 0.0)}


def test_built_in_validation_suite_passes():
    ca.run_validation_tests()


# ---------------------------------------------------------------------------
# Test integracyjny: pełne przeszukiwanie małego, syntetycznego katalogu
# ---------------------------------------------------------------------------

def _sat(norad, name, category, r, v, mass, radius):
    return {
        "norad_id": norad, "name": name, "category": category,
        "position": {"x": r[0], "y": r[1], "z": r[2]},
        "velocity": {"vx": v[0], "vy": v[1], "vz": v[2]},
        "mass_kg": mass, "radius_m": radius,
    }


@pytest.fixture
def screening(tmp_path, monkeypatch):
    # Każda para jest lustrzanym odbiciem (obrót o 180° wokół osi) -> pole grawitacyjne
    # zachowuje symetrię, więc minimalna odległość jest znana analitycznie.
    s = 750.0  # km, spotkanie po ~100 s
    catalogue = [
        # zderzenie czołowe: tor A i B przecina się dokładnie -> PHYSICAL_COLLISION
        _sat(1, "SAT-A", "active", (7000.0, -s, 0.0), (0.0, 7.5, 0.0), 500.0, 1.5),
        _sat(2, "FENGYUN 1C DEB", "debris", (7000.0, s, 0.0), (0.0, -7.5, 0.0), 1.5, 0.15),
        # przelot z odstępem ~1 km w osi z -> CRITICAL_CONJUNCTION
        _sat(3, "SAT-C", "active", (-7000.0, -s, 0.5), (0.0, 7.5, 0.0), 500.0, 1.5),
        _sat(4, "SAT-D", "active", (-7000.0, s, -0.5), (0.0, -7.5, 0.0), 500.0, 1.5),
        # para tej samej konstelacji -> pominięta
        _sat(5, "STARLINK-1", "active", (-s, 7000.0, 0.0), (7.5, 0.0, 0.0), 500.0, 1.5),
        _sat(6, "STARLINK-2", "active", (s, 7000.0, 0.0), (-7.5, 0.0, 0.0), 500.0, 1.5),
        # samotny obiekt daleko od wszystkich
        _sat(7, "LONELY", "active", (0.0, 0.0, 7200.0), (7.4, 0.0, 0.0), 500.0, 1.5),
    ]
    inp = tmp_path / "satellites_state.json"
    inp.write_text(json.dumps(catalogue))
    for name, fname in [("INPUT_FILE", inp), ("OUTPUT_JSON", tmp_path / "ca.json"),
                        ("REPORT_JSON", tmp_path / "report.json"), ("REPORT_CSV", tmp_path / "t.csv")]:
        monkeypatch.setattr(ca, name, fname)
    encounters = ca.run_catalogue_screening()
    report = json.loads((tmp_path / "report.json").read_text())
    return encounters, report


def test_screening_finds_expected_events(screening):
    encounters, _ = screening
    by_pair = {(e["norad_id_1"], e["norad_id_2"]): e for e in encounters}
    assert set(by_pair) == {(1, 2), (3, 4)}

    hit = by_pair[(1, 2)]
    assert hit["event_type"] == "PHYSICAL_COLLISION"
    assert hit["miss_distance_km"] < 1e-3
    assert hit["relative_speed_kms"] == pytest.approx(15.0, rel=0.02)
    assert hit["tca_seconds_from_epoch"] == pytest.approx(100.0, abs=2.0)
    # niezależna referencja: dokładne całkowanie toru A (scipy) do przecięcia płaszczyzny y = 0
    from scipy.integrate import solve_ivp

    def rhs(_t, s):
        r = math.sqrt(s[0] ** 2 + s[1] ** 2 + s[2] ** 2)
        k = -ca.MU_EARTH / r ** 3
        return [s[3], s[4], s[5], k * s[0], k * s[1], k * s[2]]

    cross_y = lambda _t, s: s[1]  # noqa: E731
    cross_y.terminal = True
    sol = solve_ivp(rhs, (0, 600), [7000.0, -750.0, 0.0, 0.0, 7.5, 0.0], events=cross_y,
                    rtol=1e-11, atol=1e-9)
    t_ref = sol.t_events[0][0]
    alt_ref = math.sqrt(sum(c ** 2 for c in sol.y_events[0][0][:3])) - ca.EARTH_RADIUS_KM
    assert hit["tca_seconds_from_epoch"] == pytest.approx(t_ref, abs=0.05)
    assert hit["encounter_altitude_km"] == pytest.approx(alt_ref, abs=0.5)
    assert hit["debris_involved"] is True
    assert hit["would_be_catastrophic"] is True

    near = by_pair[(3, 4)]
    assert near["event_type"] == "CRITICAL_CONJUNCTION"
    assert 0.95 < near["miss_distance_km"] < 1.0
    assert near["debris_involved"] is False


def test_screening_report_counts(screening):
    encounters, report = screening
    assert report["input_objects_count"] == 7
    assert report["unique_close_approaches_found"] == len(encounters) == 2
    assert report["physical_collisions_count"] == 1
    assert report["critical_conjunctions_under_2km"] == 1
    assert report["ignored_same_constellation_pairs"] > 0
    assert [e["miss_distance_km"] for e in encounters] == sorted(e["miss_distance_km"] for e in encounters)
