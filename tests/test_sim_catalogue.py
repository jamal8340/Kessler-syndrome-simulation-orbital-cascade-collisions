"""Engine checks on the real shared dataset (slower: loads ~18.5k objects)."""

import json

import numpy as np
import pytest

from kessler_sim import Simulation, load_scenario
from kessler_sim import orbits as ob
from kessler_sim.catalogue import DEFAULT_DATASET, DEFAULT_SUMMARY, DEFAULT_TLE_DIR, read_tle_records


@pytest.fixture(scope="module")
def sim():
    return Simulation.from_catalogue()


def test_whole_dataset_loads(sim):
    dataset = json.loads(DEFAULT_DATASET.read_text())
    assert sim.load_report.loaded == len(dataset) == len(sim.table)
    assert not sim.load_report.missing_tle and not sim.load_report.sgp4_failed
    assert sim.epoch_jd == json.loads(DEFAULT_SUMMARY.read_text())["propagation_epoch_jd"]
    assert set(sim.table.ids.tolist()) == {o["norad_id"] for o in dataset}


def _sgp4_error_km(sim, rows, dt_s):
    tles = read_tle_records(DEFAULT_TLE_DIR)
    r, _ = sim.table.states(dt_s, rows)
    jd = sim.epoch_jd + dt_s / 86400.0
    ref, ok = [], []
    for row in rows:
        err, pos, _ = tles[int(sim.table.ids[row])].sgp4(jd, 0.0)
        ok.append(err == 0)
        ref.append(pos)
    ok = np.array(ok)
    return np.linalg.norm(r[ok] - np.array(ref)[ok], axis=1)


@pytest.mark.parametrize("dt_s,median_km,p95_km", [(0.0, 15.0, 25.0), (86400.0, 15.0, 30.0),
                                                   (3 * 86400.0, 15.0, 40.0)])
def test_low_earth_orbit_accuracy_vs_sgp4(sim, dt_s, median_km, p95_km):
    """Documented accuracy: ~10 km in LEO for days (no short-period J2 terms)."""
    leo = np.flatnonzero(ob.apogee_radius(sim.table.elements) - ob.R_MEAN < 2000.0)
    rows = leo[:: max(1, len(leo) // 400)]
    err = _sgp4_error_km(sim, rows, dt_s)
    assert np.median(err) < median_km
    assert np.percentile(err, 95) < p95_km


def test_snapshot_of_whole_catalogue(sim):
    snap = sim.state(3600.0)
    assert len(snap.ids) == len(sim.table)
    alt = np.linalg.norm(snap.r_km, axis=1) - ob.R_MEAN
    assert np.all(np.isfinite(snap.r_km)) and alt.min() > 80.0


def test_closest_approach_scenario_breaks_up_a_real_pair():
    sim = load_scenario("closest_approach_breakup")
    (ev,) = sim.events
    assert ev.separation_km < 10.0
    assert ev.n_fragments > 0
    assert len(sim.state(ev.t + 60).ids) == len(sim.table) - 2
