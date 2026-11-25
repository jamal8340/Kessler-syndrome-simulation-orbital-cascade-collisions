import json
import math

import numpy as np
import pytest

from kessler_sim import BreakupParams, load_scenario
from kessler_sim import protocol
from kessler_sim.earth import earth_rotation_angle, earth_rotation_at
from test_sim_engine import T_HIT, crossing_elements, make_sim


def three_body_sim():
    return make_sim([(1, "A", "active", 900, 1.5, crossing_elements(86.0, False)),
                     (2, "B", "active", 560, 1.2, crossing_elements(74.0, True)),
                     (3, "C", "active", 500, 1.0, crossing_elements(60.0, False)),
                     (4, "D", "debris", 1.5, 0.15, crossing_elements(55.0, True))])


# -- JSON ----------------------------------------------------------------------------

def test_json_has_no_infinity_or_nan():
    text = protocol.dumps({"a": float("inf"), "b": [np.nan, -np.inf, 1.5], "c": np.float32(2.0),
                           "d": np.array([1, 2]), "e": (np.int64(3), np.bool_(True))})
    assert json.loads(text) == {"a": None, "b": [None, None, 1.5], "c": 2.0, "d": [1, 2], "e": [3, True]}
    assert "Infinity" not in text and "NaN" not in text


def test_catalogue_uses_null_for_infinite_times():
    sim = make_sim()
    entry = protocol.catalogue_entry(sim, 0)
    assert entry["t_created"] is None and entry["t_destroyed"] is None
    ev = sim.trigger_breakup(1, 2, T_HIT)
    entry = protocol.catalogue_entry(sim, 0)
    assert entry["t_destroyed"] == T_HIT
    frag = protocol.catalogue_entry(sim, sim.table.row(int(ev.fragment_ids[0])))
    assert frag["t_created"] == T_HIT and frag["t_destroyed"] is None
    assert frag["category"] == "fragment" and frag["parent_event"] == 0 and frag["parent_ids"] == [1, 2]


def test_session_and_event_messages_are_valid_json():
    sim = load_scenario("iridium_cosmos_2009")
    hello = json.loads(protocol.dumps(protocol.session(sim, "hello", 7)))
    assert hello["type"] == "hello" and hello["request_id"] == 7
    assert hello["frame"] == "TEME" and hello["epoch_utc"].startswith("2009-02-10T16:46:00")
    assert len(hello["catalogue"]) == len(sim.table)
    (ev,) = hello["events"]
    assert ev["type"] == "breakup_event" and ev["remnant_id"] is None
    n = ev["n_fragments"]
    assert len(ev["fragment_ids"]) == len(ev["fragment_parent_ids"]) == len(ev["fragment_dv_kms"]) == n
    assert set(ev["fragment_parent_ids"]) == {24946, 22675}


# -- binary state ------------------------------------------------------------------

def test_state_round_trip_with_nan_for_missing_objects():
    sim = make_sim()
    sim.trigger_breakup(1, 2, T_HIT)
    for t, parents_exist in ((T_HIT - 1, True), (T_HIT + 1, False)):
        state = protocol.unpack_state(protocol.state_message(sim, t))
        assert state["t"] == t and state["catalogue_version"] == sim.table.version
        assert state["earth_rotation_rad"] == pytest.approx(earth_rotation_at(sim, t))
        pos = state["positions"]
        assert pos.shape == (len(sim.table), 3)
        exists = ~np.isnan(pos[:, 0])
        assert exists[:2].tolist() == [parents_exist] * 2
        assert exists[2:].all() != parents_exist
        expected = sim.state(t)
        assert np.allclose(pos[exists], expected.r_km, atol=1e-2)   # float32
        assert sim.table.ids[exists].tolist() == expected.ids.tolist()


def test_state_header_layout():
    data = protocol.pack_state(12.5, 1.0, 3, np.zeros((2, 3), np.float32))
    assert len(data) == 8 + 8 + 4 + 4 + 2 * 3 * 4


# -- Earth rotation ----------------------------------------------------------------

def test_earth_rotation_reference_value():
    # Vallado, Fundamentals of Astrodynamics, example 3-5: 1992-08-20 12:14 UT1 -> GMST 152.578788 deg
    assert math.degrees(earth_rotation_angle(2448855.009722)) == pytest.approx(152.578788, abs=1e-4)


def test_earth_rotation_one_sidereal_day():
    jd = 2461357.45
    a0, a1 = earth_rotation_angle(jd), earth_rotation_angle(jd + 0.99726957)
    assert abs((a1 - a0 + math.pi) % (2 * math.pi) - math.pi) < 1e-5


# -- change tracking: applying updates like Unity reproduces the engine ---------------

class FakeUnity:
    """Applies catalogue updates exactly as TRANSPORT_ARCHITECTURE.md prescribes."""

    def __init__(self, sim):
        self.entries = protocol.catalogue(sim)
        self.events = {e.event_id for e in sim.events}

    def apply(self, changes, sim):
        update = json.loads(protocol.dumps(protocol.catalogue_update(sim, changes)))
        removed = set(update["removed"])
        self.entries = [e for e in self.entries if e["id"] not in removed]       # 1. remove
        by_id = {e["id"]: e for e in update["updated"]}
        self.entries = [by_id.get(e["id"], e) for e in self.entries]             # 2. update
        self.entries += update["added"]                                          # 3. append
        self.events -= set(changes.removed_events)
        self.events |= set(changes.added_events)
        assert len(self.entries) == update["count"]


def check_in_sync(unity, sim):
    expected = json.loads(protocol.dumps(protocol.catalogue(sim)))
    assert unity.entries == expected
    assert unity.events == {e.event_id for e in sim.events}


def test_breakup_changes():
    sim = make_sim()
    sim.take_changes()
    ev = sim.trigger_breakup(1, 2, T_HIT)
    changes = sim.take_changes()
    assert changes.removed == [] and changes.updated == [1, 2]
    assert changes.added == ev.created_ids and changes.added_events == [0]
    assert not sim.take_changes()


def test_rewrite_reuses_ids_and_updates_apply_cleanly():
    sim = three_body_sim()
    unity = FakeUnity(sim)
    sim.take_changes()

    first = sim.trigger_breakup(1, 2, T_HIT)
    unity.apply(sim.take_changes(), sim)
    check_in_sync(unity, sim)

    second = sim.trigger_breakup(1, 3, T_HIT - 50)            # earlier: erases the first breakup
    changes = sim.take_changes()
    assert changes.removed_events == [0] and changes.added_events == [0]
    assert set(changes.removed) == set(first.created_ids)
    assert set(changes.added) == set(second.created_ids)
    assert set(changes.removed) & set(changes.added)          # same IDs removed and re-added
    assert 2 in changes.updated                                # B was revived
    unity.apply(changes, sim)
    check_in_sync(unity, sim)


def test_random_sequence_of_breakups_and_rewrites_stays_in_sync():
    rng = np.random.default_rng(0)
    sim = three_body_sim()
    unity = FakeUnity(sim)
    sim.take_changes()
    for _ in range(25):
        t = float(rng.uniform(T_HIT - 200, T_HIT + 200))
        alive = sim.state(t).ids
        if len(alive) < 2:
            continue
        a, b = rng.choice(alive, 2, replace=False)
        sim.trigger_breakup(int(a), int(b), t, BreakupParams(max_fragments=30))
        if rng.random() < 0.5:                                # sometimes batch two operations
            continue
        unity.apply(sim.take_changes(), sim)
        check_in_sync(unity, sim)
    unity.apply(sim.take_changes(), sim)
    check_in_sync(unity, sim)
