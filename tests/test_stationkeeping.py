"""Ideal compensation for intact active satellites; passive debris is unchanged."""
import numpy as np
import pytest

from kessler_sim import DragConfig, Simulation, BreakupParams, drag, orbits

DAY = 86400.
YEAR = 365.25 * DAY


def mixed():
    sim = Simulation()
    n = 6
    z = np.zeros(n)
    a = orbits.R_MEAN + np.array([400., 400., 400., 400., 2500., 149.])
    sim.table.add(range(1, n+1), ['active', 'debris', 'fragment', 'rocket', 'high', 'low'],
                  ['active', 'debris', 'fragment', 'rocket_body', 'active', 'active'],
                  100., 1., orbits.Elements(a, z, z+.7, z, z, z, z+1e-12),
                  0., am_m2_per_kg=.02)
    return sim


def test_maintained_vs_passive_and_reentry_threshold():
    sim = mixed()
    sim.enable_drag()
    sim.advance_drag(50*YEAR)
    assert sim.state(50*YEAR).ids.tolist() == [1, 5]
    assert set(sim.reentries) == {2, 3, 4, 6}
    el = sim.table.elements_at(50*YEAR, [0, 4])
    np.testing.assert_array_equal(el.a, sim.table.elements.a[[0, 4]])
    assert np.all(el.ndot == 0)
    assert sim.table.orbit_history.latest[0, drag.END] == np.inf


def test_maintained_phase_j2_and_rewind():
    sim = mixed()
    t0 = 12345.
    initial = sim.table.elements_at(t0, [0])
    old = sim.table.states(t0-1, [0])
    start = sim.table.states(t0, [0])
    sim.enable_drag(t=t0)
    sim.advance_drag(YEAR)
    initial.ndot[:] = 0
    expected = orbits.elements_to_state(initial, YEAR-t0)
    for got, want in zip(sim.table.states(YEAR, [0]), expected):
        np.testing.assert_allclose(got, want, atol=1e-8)
    for t, original in [(t0-1, old), (t0, start)]:
        for got, want in zip(sim.table.states(t, [0]), original):
            np.testing.assert_allclose(got, want, atol=1e-8)
    sim.truncate_after(t0-1)
    assert sim.drag_config is None
    assert sim.table.orbit_history.segment_count == 0
    assert sim.table.elements_at(YEAR, [0]).ndot[0] == 1e-12


def test_passive_switch_restores_original_drag():
    sim = mixed()
    sim.enable_drag(DragConfig(stationkeeping=False))
    sim.advance_drag(YEAR)
    assert sim.reentries[1] == pytest.approx(sim.reentries[2])
    assert sim.table.elements_at(YEAR, [4]).ndot[0] == 1e-12


def test_active_remnant_and_fragments_are_passive_after_collision():
    sim = Simulation()
    z = np.zeros(2)
    sim.table.add([1,2], ['satellite','bolt'], ['active','debris'], [1000., .02], [1.5,.05],
                  orbits.Elements(z+orbits.R_MEAN+550., z, np.radians([86.,74.]), z, z, z), 0.)
    sim.enable_drag()
    sim.advance_drag(DAY)
    event = sim.trigger_breakup(1,2,DAY, BreakupParams(impact_speed_kms=10.))
    assert not event.catastrophic and event.survivor_id is not None
    row = sim.table.row(event.survivor_id)
    assert sim.table.category[row] == 'active'  # legacy UI category retained
    a0 = sim.table.elements_at(DAY,[row]).a[0]
    sim.advance_drag(2*DAY)
    assert sim.table.is_alive(row, 2*DAY)
    assert sim.table.elements_at(2*DAY,[row]).a[0] < a0
    assert np.isfinite(sim.table.orbit_history.latest[row,drag.END])
    for oid in event.fragment_ids:
        assert not np.isposinf(sim.table.orbit_history.latest[sim.table.row(oid),drag.END])
    sim.truncate_after(DAY-1)
    assert sim.table.is_alive(sim.table.row(1),DAY-1)
    assert not sim.table.has(event.survivor_id)
    sim.advance_drag(2*DAY)
    assert sim.table.elements_at(2*DAY,[0]).a[0] == sim.table.elements.a[0]
