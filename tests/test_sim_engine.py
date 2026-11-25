import numpy as np
import pytest

from kessler_sim import BreakupConfig, BreakupParams, Simulation, available_scenarios, load_scenario
from kessler_sim import orbits as ob
from kessler_sim.analysis import cloud_stats
from kessler_sim.clock import Clock
from kessler_sim.state import ObjectTable

T_HIT = 300.0
ALT = 789.0


def crossing_elements(inc_deg, descending, lat_deg=50.0, lon_deg=0.0, radius=None):
    """Circular orbit passing through the inertial point (lat, lon, radius) at its t_ref."""
    a = radius if radius is not None else ob.R_MEAN + ALT
    i, lat = np.radians(inc_deg), np.radians(lat_deg)
    u = np.arcsin(np.sin(lat) / np.sin(i))
    if descending:
        u = np.pi - u
    raan = np.radians(lon_deg) - np.arctan2(np.cos(i) * np.sin(u), np.cos(u))
    return ob.Elements(np.array([a]), np.array([0.0]), np.array([i]), np.array([raan % (2 * np.pi)]),
                       np.array([0.0]), np.array([u]))


def make_sim(objects=None, **config):
    """Objects: list of (id, name, category, mass, radius, Elements). Default: one crossing pair."""
    if objects is None:
        objects = [
            (1, "SAT-A", "active", 900.0, 1.5, crossing_elements(86.0, False)),
            (2, "SAT-B", "active", 560.0, 1.2, crossing_elements(74.0, True)),
        ]
    table = ObjectTable()
    for oid, name, cat, mass, radius, el, *t_ref in objects:
        table.add([oid], name, cat, mass, radius, el, t_ref=t_ref[0] if t_ref else T_HIT)
    return Simulation(table, epoch_jd=2461357.5, config=BreakupConfig(**config))


def momentum_check(sim, ev):
    """Momentum of fragments (+ unassigned mass moving with its parent) equals the parents' momentum."""
    tb = sim.table
    (va, vb) = tb.states(ev.t, [tb.row(ev.id_a), tb.row(ev.id_b)])[1]
    rows = [tb.row(i) for i in ev.fragment_ids]
    _, v = tb.states(ev.t, rows)
    m = tb.mass_kg[rows]
    p = (m[:, None] * v).sum(axis=0)
    for pid, vp, mp in ((ev.id_a, va, ev.mass_a_kg), (ev.id_b, vb, ev.mass_b_kg)):
        assigned = m[ev.fragment_origin_ids == pid].sum()
        p += (mp - assigned) * vp
    return p, ev.mass_a_kg * va + ev.mass_b_kg * vb


# -- clock -------------------------------------------------------------------------

def test_clock_play_pause_speed_seek():
    c = Clock(speed=60.0)
    assert c.advance(1.0) == 0.0            # paused by default
    c.play()
    assert c.advance(1.0) == 60.0
    c.set_speed(-120.0)
    assert c.advance(0.25) == 30.0
    assert c.advance(10.0) == 0.0 and not c.running   # clamped at epoch, stops
    assert c.seek(-5) == 0.0
    assert c.seek(1234.5) == 1234.5
    c.pause()
    assert c.advance(5.0) == 1234.5


# -- object table ------------------------------------------------------------------

def test_table_rejects_duplicate_ids():
    sim = make_sim()
    with pytest.raises(ValueError):
        sim.table.add([1], "X", "active", 1.0, 1.0, crossing_elements(50.0, False), t_ref=0.0)


def test_state_returns_alive_objects():
    sim = make_sim()
    snap = sim.state(T_HIT)
    assert sorted(snap.ids.tolist()) == [1, 2]
    assert np.allclose(snap.r_km[0], snap.r_km[1], atol=1e-6)     # they meet at T_HIT
    assert snap.r_km.shape == snap.v_kms.shape == (2, 3)


def test_constellation_tag():
    els = [crossing_elements(53.0, False)] * 3
    sim = make_sim([(10, "STARLINK-1", "active", 300, 1, els[0]), (11, "STARLINK DEB", "debris", 1, 0.1, els[1]),
                    (12, "ISS (ZARYA)", "active", 4e5, 50, els[2])])
    assert sim.table.constellation.tolist() == ["STARLINK", "", ""]


# -- breakups ----------------------------------------------------------------------

def test_catastrophic_breakup_creates_two_clouds():
    sim = make_sim()
    ev = sim.trigger_breakup(1, 2, T_HIT)
    assert ev.catastrophic and ev.survivor_id is None
    assert ev.separation_km < 1e-6
    assert ev.n_fragments > 500
    tb = sim.table
    for pid, inc in ((1, 86.0), (2, 74.0)):
        rows = [tb.row(i) for i in ev.fragment_ids[ev.fragment_origin_ids == pid]]
        assert len(rows) > 100
        assert abs(np.degrees(tb.elements.i[rows]).mean() - inc) < 0.5


def test_parents_exist_until_breakup_and_fragments_after():
    sim = make_sim()
    ev = sim.trigger_breakup(1, 2, T_HIT)
    before, at, after = sim.state(T_HIT - 1), sim.state(T_HIT), sim.state(T_HIT + 3600)
    assert sorted(before.ids.tolist()) == [1, 2]
    for snap in (at, after):
        assert 1 not in snap.ids and 2 not in snap.ids
        assert set(snap.ids.tolist()) == set(ev.fragment_ids.tolist())


def test_fragments_start_at_collision_point_with_parent_velocity_plus_kick():
    sim = make_sim()
    (va, vb) = sim.table.states(T_HIT, [0, 1])[1]
    ev = sim.trigger_breakup(1, 2, T_HIT)
    rows = [sim.table.row(i) for i in ev.fragment_ids]
    r, v = sim.table.states(T_HIT, rows)
    assert np.allclose(r, ev.position_km, atol=1e-5)
    parent_v = np.where((ev.fragment_origin_ids == 1)[:, None], va, vb)
    assert np.allclose(v, parent_v + ev.fragment_dv_kms, atol=1e-8)


def test_momentum_conserved_per_parent():
    # keep every fragment: no re-entry filter, and l_min large enough that none escapes
    sim = make_sim(reentry_altitude_km=-6000.0)
    ev = sim.trigger_breakup(1, 2, T_HIT, BreakupParams(l_min_m=0.3))
    assert ev.reentered_count == 0 and ev.escaped_count == 0
    p_after, p_before = momentum_check(sim, ev)
    assert np.allclose(p_after, p_before, rtol=1e-9, atol=1e-6)


def test_mass_accounted_for():
    sim = make_sim()
    ev = sim.trigger_breakup(1, 2, T_HIT)
    total = ev.tracked_mass_kg + ev.untracked_mass_kg + ev.lost_mass_kg
    assert total == pytest.approx(900.0 + 560.0, rel=1e-9)
    assert sim.table.mass_kg[[sim.table.row(i) for i in ev.fragment_ids]].sum() == pytest.approx(ev.tracked_mass_kg)


def test_reentering_fragments_are_not_added():
    sim = make_sim(reentry_altitude_km=600.0)
    ev = sim.trigger_breakup(1, 2, T_HIT)
    assert ev.reentered_count > 0
    rows = [sim.table.row(i) for i in ev.fragment_ids]
    perigee = ob.perigee_radius(sim.table.elements.take(rows)) - ob.R_MEAN
    assert perigee.min() > 600.0


def test_non_catastrophic_breakup_leaves_remnant():
    small = crossing_elements(74.0, True)
    sim = make_sim([(1, "BIG", "active", 1000.0, 1.5, crossing_elements(86.0, False)),
                    (2, "BOLT", "debris", 0.02, 0.05, small)])
    ev = sim.trigger_breakup(1, 2, T_HIT)
    assert not ev.catastrophic
    assert ev.survivor_id is not None
    info = sim.object_info(ev.survivor_id, T_HIT + 1)
    assert info["alive"] and info["name"] == "BIG [REMNANT]"
    assert info["mass_kg"] + ev.tracked_mass_kg + ev.lost_mass_kg == pytest.approx(1000.02)
    assert abs(info["orbit"]["i_deg"] - 86.0) < 0.5            # follows the target's orbit


def test_breakup_errors():
    sim = make_sim()
    with pytest.raises(ValueError):
        sim.trigger_breakup(1, 1, T_HIT)
    with pytest.raises(KeyError):
        sim.trigger_breakup(1, 99, T_HIT)
    sim.trigger_breakup(1, 2, T_HIT)
    frag = int(sim.events[0].fragment_ids[0])
    with pytest.raises(ValueError):                 # parent no longer exists after the breakup
        sim.trigger_breakup(1, frag, T_HIT + 10)
    with pytest.raises(ValueError):                 # fragment does not exist yet
        sim.trigger_breakup(frag, int(sim.events[0].fragment_ids[1]), T_HIT - 10)


def test_overrides_impact_speed_masses_and_cap():
    sim = make_sim()
    ev = sim.trigger_breakup(1, 2, T_HIT, BreakupParams(impact_speed_kms=3.0, mass_a_kg=50.0, mass_b_kg=40.0,
                                                         max_fragments=100))
    assert ev.v_rel_kms == pytest.approx(3.0)
    assert (ev.mass_a_kg, ev.mass_b_kg) == (50.0, 40.0)
    assert ev.n_expected <= 101 and ev.l_min_m > 0.1


def test_fragment_cap_raises_minimum_size():
    sim = make_sim(max_fragments=200)
    ev = sim.trigger_breakup(1, 2, T_HIT)
    assert ev.n_expected <= 201
    assert ev.l_min_m > 0.1
    assert ev.fragment_lc_m.min() >= ev.l_min_m - 1e-12


def test_breakup_reproducible_with_seed():
    a, b = make_sim(), make_sim()
    ea, eb = a.trigger_breakup(1, 2, T_HIT), b.trigger_breakup(1, 2, T_HIT)
    assert np.array_equal(ea.fragment_ids, eb.fragment_ids)
    assert np.array_equal(a.table.elements.a, b.table.elements.a)
    c = make_sim(seed=7)
    assert not np.array_equal(c.trigger_breakup(1, 2, T_HIT).fragment_mass_kg, ea.fragment_mass_kg)


def test_new_breakup_in_the_past_rewrites_history():
    els = [crossing_elements(86.0, False), crossing_elements(74.0, True), crossing_elements(30.0, False, 20.0)]
    sim = make_sim([(1, "A", "active", 900, 1.5, els[0]), (2, "B", "active", 560, 1.2, els[1]),
                    (3, "C", "active", 500, 1.0, els[2])])
    first = sim.trigger_breakup(1, 2, T_HIT)
    n_first = len(sim.table)
    second = sim.trigger_breakup(1, 3, T_HIT - 100)   # earlier: object 1 hits object 3 instead
    assert len(sim.events) == 1 and sim.events[0].id_b == 3
    # fragment IDs are reused (deterministic), but every fragment now belongs to the new event
    frag_rows = sim.fragment_rows(0)
    assert set(sim.table.ids[frag_rows].tolist()) == set(second.fragment_ids.tolist())
    assert all(sim.table.parent_ids[r] == (1, 3) for r in frag_rows)
    assert sim.table.is_alive(sim.table.row(2), T_HIT + 10)   # B was never hit in the new history
    assert len(sim.table) == 3 + second.n_fragments != n_first
    assert first.fragment_ids[0] == second.fragment_ids[0]


def test_rewind_then_same_breakup_gives_same_fragments():
    sim = make_sim()
    first = sim.trigger_breakup(1, 2, T_HIT)
    ids, masses = first.fragment_ids.copy(), first.fragment_mass_kg.copy()
    sim.truncate_after(T_HIT - 1)
    assert sim.events == [] and len(sim.table) == 2
    again = sim.trigger_breakup(1, 2, T_HIT)
    assert np.array_equal(again.fragment_ids, ids) and np.array_equal(again.fragment_mass_kg, masses)


def test_object_info_fields():
    sim = make_sim()
    info = sim.object_info(1, 0.0)
    assert info["alive"] and info["parent_event"] is None and info["t_created"] is None
    assert info["altitude_km"] == pytest.approx(ALT, abs=1e-6)
    assert info["orbit"]["i_deg"] == pytest.approx(86.0)
    assert info["orbit"]["period_min"] == pytest.approx(100.4, abs=0.2)


def test_object_info_has_no_position_when_object_does_not_exist():
    sim = make_sim()
    ev = sim.trigger_breakup(1, 2, T_HIT)
    frag = int(ev.fragment_ids[0])
    for oid, t in ((1, T_HIT + 100), (1, T_HIT), (frag, 0.0), (frag, T_HIT - 0.001)):
        info = sim.object_info(oid, t)
        assert not info["alive"]
        assert info["position_km"] is None and info["velocity_kms"] is None and info["altitude_km"] is None
    assert sim.object_info(1, T_HIT - 1)["position_km"] is not None        # parent before the breakup
    assert sim.object_info(frag, T_HIT)["position_km"] is not None         # fragment from its creation
    assert sim.object_info(frag, 0.0)["t_created"] == T_HIT
    assert sim.object_info(1, T_HIT + 100)["t_destroyed"] == T_HIT


# -- close approaches --------------------------------------------------------------

def test_screening_finds_crossing_pair():
    sim = make_sim()
    found = sim.find_close_approaches(0.0, 600.0)
    assert len(found) == 1
    hit = found[0]
    assert (hit["id_1"], hit["id_2"]) == (1, 2)
    assert hit["tca_s"] == pytest.approx(T_HIT, abs=0.01)
    assert hit["miss_distance_km"] < 1e-3
    assert hit["event_type"] == "PHYSICAL_COLLISION"
    assert hit["altitude_km"] == pytest.approx(ALT, abs=0.1)
    assert hit["would_be_catastrophic"]


def test_screening_skips_same_constellation_and_dead_objects():
    sim = make_sim([(1, "STARLINK-1", "active", 300, 1.5, crossing_elements(86.0, False)),
                    (2, "STARLINK-2", "active", 300, 1.5, crossing_elements(74.0, True))])
    assert sim.find_close_approaches(0.0, 600.0) == []
    assert len(sim.find_close_approaches(0.0, 600.0, skip_same_constellation=False)) == 1

    sim = make_sim()
    sim.trigger_breakup(1, 2, T_HIT - 100)            # destroyed before they would meet
    found = sim.find_close_approaches(0.0, 600.0)
    assert not any({e["id_1"], e["id_2"]} == {1, 2} for e in found)


def test_screening_ignores_sibling_fragments_but_sees_other_objects():
    # C crosses the point where parent A (and so A's debris cloud) is 30 s after the collision
    a_el = crossing_elements(86.0, False)
    (p,), _ = ob.elements_to_state(a_el, 30.0)
    lat, lon = np.degrees(np.arcsin(p[2] / np.linalg.norm(p))), np.degrees(np.arctan2(p[1], p[0]))
    third = crossing_elements(60.0, True, lat, lon, np.linalg.norm(p))
    sim = make_sim([(1, "A", "active", 900, 1.5, a_el),
                    (2, "B", "active", 560, 1.2, crossing_elements(74.0, True)),
                    (3, "C", "active", 500, 1.0, third, T_HIT + 30)])
    ev = sim.trigger_breakup(1, 2, T_HIT)
    found = sim.find_close_approaches(T_HIT, 120.0)
    frag_ids = set(ev.fragment_ids.tolist())
    assert found, "object C should meet the debris cloud"
    for e in found:
        assert not ({e["id_1"], e["id_2"]} <= frag_ids)      # never fragment vs sibling fragment
        assert 3 in (e["id_1"], e["id_2"])


# -- scenarios and analysis --------------------------------------------------------

def test_available_scenarios():
    assert {"catalogue", "closest_approach_breakup", "iridium_cosmos_2009"} <= set(available_scenarios())
    with pytest.raises(FileNotFoundError):
        load_scenario("does_not_exist")


def test_iridium_cosmos_scenario():
    sim = load_scenario("iridium_cosmos_2009")
    assert sim.epoch_utc.strftime("%Y-%m-%d %H:%M") == "2009-02-10 16:46"
    (ev,) = sim.events
    assert ev.t == 600.0 and ev.separation_km < 1e-6
    assert ev.v_rel_kms == pytest.approx(11.6, abs=0.1)
    assert ev.catastrophic
    assert 900 < ev.n_fragments < 1300
    counts = ev.summary()["fragments_per_parent"]
    assert counts[24946] > 300 and counts[22675] > 300


def test_cloud_spreads_along_the_orbit():
    sim = load_scenario("iridium_cosmos_2009")
    spread = [cloud_stats(sim, 0, 600.0 + dt)["arc_spread_p95_deg"] for dt in (0, 600, 6 * 3600, 86400)]
    assert spread[0] < 1e-3
    assert spread == sorted(spread)
    assert spread[-1] > 90.0
    stats = cloud_stats(sim, 0, 600.0 + 30 * 86400)
    assert stats["per_parent"][22675]["raan_std_deg"] > 2.0      # J2 spreads the orbital planes
