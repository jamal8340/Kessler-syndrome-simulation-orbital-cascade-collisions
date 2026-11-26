"""Simulation engine: the single entry point the server (and scripts, tests) use.

It owns the object table, the clock and the event log. Positions are computed on
demand for any time t, so play/pause/seek are trivial and runs are reproducible:
the state is fully determined by the dataset, the seed and the list of breakups.

Rewriting history: triggering a breakup at time t discards every breakup that
happened after t (their fragments are removed and their parents come back).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import drag, orbits, screening
from .breakup import BreakupConfig, BreakupEvent, BreakupParams, apply_breakup
from .catalogue import jd_to_datetime, load_catalogue
from .clock import Clock
from .state import ObjectTable


@dataclass
class Snapshot:
    t: float
    version: int
    ids: np.ndarray        # (n,)
    r_km: np.ndarray       # (n, 3)
    v_kms: np.ndarray      # (n, 3)


@dataclass
class CatalogueChanges:
    """What changed in the object table since the last ``take_changes()``.

    Apply in this order: remove ``removed``, update the remaining ``updated`` objects
    (their lifespan changed), then append ``added`` in order. IDs can be reused after a
    rewrite of history, so the same ID may appear in ``removed`` and ``added``.
    """

    removed: list[int] = field(default_factory=list)
    updated: list[int] = field(default_factory=list)
    added: list[int] = field(default_factory=list)
    removed_events: list[int] = field(default_factory=list)
    added_events: list[int] = field(default_factory=list)

    def __bool__(self):
        return bool(self.removed or self.updated or self.added or self.removed_events or self.added_events)


class Simulation:
    def __init__(self, table: ObjectTable | None = None, epoch_jd: float = 2451545.0,
                 config: BreakupConfig | None = None, name: str = ""):
        self.table = table if table is not None else ObjectTable()
        self.epoch_jd = float(epoch_jd)
        self.config = config or BreakupConfig()
        self.name = name
        self.clock = Clock()
        self.events: list[BreakupEvent] = []
        self.drag_config: drag.DragConfig | None = None
        self.drag_start: float | None = None
        self.drag_time: float | None = None
        self.reentries: dict[int, float] = {}  # drag re-entries only, not immediate breakup losses
        self._removed: list[int] = []
        self._added: list[int] = []
        self._updated: set[int] = set()
        self._removed_events: list[int] = []
        self._added_events: list[int] = []

    @classmethod
    def from_catalogue(cls, config: BreakupConfig | None = None, name="catalogue", **load_kwargs):
        table, epoch_jd, report = load_catalogue(**load_kwargs)
        sim = cls(table, epoch_jd, config, name)
        sim.load_report = report
        return sim

    # -- time ------------------------------------------------------------------

    @property
    def epoch_utc(self):
        return jd_to_datetime(self.epoch_jd)

    @property
    def t(self) -> float:
        return self.clock.t

    def utc_at(self, t: float):
        return jd_to_datetime(self.epoch_jd + t / 86400.0)

    # -- state -------------------------------------------------------------------

    def state(self, t: float | None = None) -> Snapshot:
        """All objects alive at ``t`` (default: clock time) with their positions."""
        t = self.clock.t if t is None else float(t)
        rows = self.table.alive_rows(t)
        r, v = self.table.states(t, rows)
        return Snapshot(t, self.table.version, self.table.ids[rows].copy(), r, v)

    def object_info(self, object_id: int, t: float | None = None) -> dict:
        t = self.clock.t if t is None else float(t)
        tb = self.table
        row = tb.row(object_id)
        el = tb.elements_at(t, [row])
        (r,), (v,) = tb.states(t, [row])
        alive = tb.is_alive(row, t)
        info = {
            "id": int(tb.ids[row]), "name": str(tb.names[row]), "category": str(tb.category[row]),
            "mass_kg": float(tb.mass_kg[row]), "radius_m": float(tb.radius_m[row]),
            "am_m2_per_kg": None if np.isnan(tb.am_m2_per_kg[row]) else float(tb.am_m2_per_kg[row]),
            "alive": alive,
            "t_created": None if np.isinf(tb.t_created[row]) else float(tb.t_created[row]),
            "t_destroyed": None if np.isinf(tb.t_destroyed[row]) else float(tb.t_destroyed[row]),
            "parent_event": None if tb.parent_event[row] < 0 else int(tb.parent_event[row]),
            "parent_ids": list(tb.parent_ids[row]),
            # no position for an object that does not exist at t (before creation / after destruction)
            "position_km": r.tolist() if alive else None,
            "velocity_kms": v.tolist() if alive else None,
            "altitude_km": float(np.linalg.norm(r) - orbits.R_MEAN) if alive else None,
            "orbit": {
                "a_km": float(el.a[0]), "e": float(el.e[0]), "i_deg": float(np.degrees(el.i[0])),
                "perigee_alt_km": float(orbits.perigee_radius(el)[0] - orbits.R_MEAN),
                "apogee_alt_km": float(orbits.apogee_radius(el)[0] - orbits.R_MEAN),
                "period_min": float(orbits.period(el)[0] / 60.0),
            },
        }
        return info

    # -- change tracking (for the server) ----------------------------------------

    def _record_removed(self, ids):
        for i in ids:
            if i in self._added:
                self._added.remove(i)       # never announced: just forget it
            else:
                self._removed.append(i)
            self._updated.discard(i)

    def take_changes(self) -> CatalogueChanges:
        added = set(self._added)
        removed_only = set(self._removed) - added
        changes = CatalogueChanges(
            removed=list(self._removed),
            updated=[i for i in sorted(self._updated) if i not in added and i not in removed_only],
            added=list(self._added),
            removed_events=list(self._removed_events),
            added_events=list(self._added_events),
        )
        self._removed, self._added, self._updated = [], [], set()
        self._removed_events, self._added_events = [], []
        return changes

    # -- optional long-term drag -------------------------------------------------

    def enable_drag(self, config: drag.DragConfig | None = None, *, t: float = 0.) -> None:
        """Enable drag after any preloaded breakups; existing scenarios remain unchanged.

        The batch runner must advance drag before sampling collisions at a later time.
        Seeking is read-only; it never runs physics or modifies the recorded history.
        """
        t = float(t)
        if not np.isfinite(t) or (self.events and t < self.events[-1].t):
            raise ValueError("drag must start at or after the last preloaded breakup")
        if self.drag_config is not None:
            raise ValueError("drag is already enabled")
        config = config or drag.DragConfig(reentry_altitude_km=self.config.reentry_altitude_km)
        if config.reentry_altitude_km != self.config.reentry_altitude_km:
            raise ValueError("drag and breakup re-entry thresholds must match")
        # Validate physical inputs before mutating the simulation.
        rows = self.table.alive_rows(t)
        drag.area_to_mass(self.table.am_m2_per_kg[rows], self.table.radius_m[rows], self.table.mass_kg[rows])
        self.drag_config, self.drag_start, self.drag_time = config, t, t
        self.advance_drag(t)

    def _record_reentries(self, rows, times):
        tb = self.table
        times = np.broadcast_to(times, (len(rows),))
        if len(rows):
            tb.destroy(rows, times)
            for row, time in zip(rows, times):
                oid = int(tb.ids[row])
                self.reentries[oid] = float(time)
                self._updated.add(oid)

    def advance_drag(self, t: float) -> list[tuple[int, float]]:
        """Generate/commit drag history through t, independently of collision sampling cadence.

        Adaptive segments may predict beyond t; predicted re-entry is not committed
        until its time is reached. Return newly committed (ID, time) pairs in time order.
        """
        t = float(t)
        cfg = self.drag_config
        if cfg is None or self.drag_time is None:
            raise ValueError("drag is not enabled")
        if not np.isfinite(t) or t < self.drag_time:
            raise ValueError("drag time must advance monotonically; truncate history first")
        tb = self.table
        before = set(self.reentries)
        threshold = orbits.R_MEAN + cfg.reentry_altitude_km
        rows = np.flatnonzero((tb.t_created <= t) & (tb.t_destroyed > self.drag_time))
        # Include scheduled births inside this interval, and new collision fragments.
        new = rows[~np.isfinite(tb.orbit_history.latest[rows, drag.START])]
        starts = np.maximum(tb.t_created[new], self.drag_time)
        for start in np.unique(starts):
            born = new[starts == start]
            el = tb.elements_at(start, born)
            peri = orbits.perigee_radius(el)
            valid = np.isfinite(el.a) & (el.a > 0) & (el.e >= 0) & (el.e < 1)
            if not valid.all():
                raise ValueError("drag requires bound, valid mean orbits")
            reentry = peri <= threshold
            self._record_reentries(born[reentry], start)
            # Only intact active objects are assumed to maintain their orbit.
            # A non-catastrophic remnant may retain category=active, but is passive.
            maintained = (cfg.stationkeeping & ~reentry &
                          (tb.category[born] == "active") & (tb.parent_event[born] == -1))
            if maintained.any():
                tb.orbit_history.append(tb.ids, born[maintained],
                                        drag.maintained_segments(el.take(maintained), start))
            use = ~reentry & ~maintained & (peri <= orbits.R_MEAN + cfg.max_perigee_altitude_km)
            if use.any():
                selected = born[use]
                am = drag.area_to_mass(tb.am_m2_per_kg[selected], tb.radius_m[selected], tb.mass_kg[selected])
                # Rebase at the start and replace catalogue ndot, never apply both models.
                segments = drag.plan_segments(el.take(use), start, am, cfg)
                tb.orbit_history.append(tb.ids, selected, segments)
        # Each iteration advances only objects whose own segment has expired.
        # Other objects retain their current segment even if collision ticks are frequent.
        while True:
            latest = tb.orbit_history.latest
            due = rows[(latest[rows, drag.END] <= t) & (tb.t_destroyed[rows] > self.drag_time)]
            if not len(due):
                break
            ends = latest[due, drag.END].copy()
            el = drag.segment_elements(latest[due], ends)
            reentry = orbits.perigee_radius(el) <= threshold + 1e-8
            self._record_reentries(due[reentry], ends[reentry])
            # Remove re-entered rows from the work list to avoid processing them again.
            rows = rows[~np.isin(rows, due[reentry])] if reentry.any() else rows
            keep = ~reentry
            if keep.any():
                selected = due[keep]
                am = drag.area_to_mass(tb.am_m2_per_kg[selected], tb.radius_m[selected], tb.mass_kg[selected])
                segments = drag.plan_segments(el.take(keep), ends[keep], am, cfg)
                tb.orbit_history.append(tb.ids, selected, segments)
        self.drag_time = t
        return sorted(((oid, time) for oid, time in self.reentries.items() if oid not in before),
                      key=lambda item: (item[1], item[0]))

    # -- events ------------------------------------------------------------------

    def truncate_after(self, t: float) -> list[BreakupEvent]:
        """Discard future drag/re-entries and breakups; restore their parents."""
        t = float(t)
        if not np.isfinite(t):
            raise ValueError("history time must be finite")
        if self.drag_time is not None and t < self.drag_time:
            self.table.orbit_history.truncate_after(self.table.ids, t)
            for oid, time in list(self.reentries.items()):
                if time > t:
                    row = self.table.row(oid)
                    self.table.revive([row])
                    self._updated.add(oid)
                    del self.reentries[oid]
            if t < self.drag_start:
                self.drag_config = self.drag_start = self.drag_time = None
            else:
                self.drag_time = t
        removed = []
        while self.events and self.events[-1].t > t:
            ev = self.events.pop()
            created_ids = [i for i in ev.created_ids if self.table.has(i)]
            parents = [self.table.row(ev.id_a), self.table.row(ev.id_b)]
            self.table.revive(parents)
            self._updated.update((ev.id_a, ev.id_b))
            if created_ids:
                self.table.remove_rows([self.table.row(i) for i in created_ids])
                self._record_removed(created_ids)
            if ev.event_id in self._added_events:
                self._added_events.remove(ev.event_id)
            else:
                self._removed_events.append(ev.event_id)
            removed.append(ev)
        return removed

    def trigger_breakup(self, id_a: int, id_b: int, t: float | None = None,
                        params: BreakupParams | None = None) -> BreakupEvent:
        t = self.clock.t if t is None else float(t)
        if not np.isfinite(t):
            raise ValueError("breakup time must be finite")
        if self.drag_time is not None and t > self.drag_time:
            raise ValueError("advance drag to the breakup time before triggering a collision")
        # validate before rewriting history, so an invalid request never discards events
        if int(id_a) == int(id_b):
            raise ValueError("an object cannot collide with itself")
        for oid in (id_a, id_b):
            if not self.table.is_alive(self.table.row(oid), t):
                raise ValueError(f"object {oid} does not exist at t = {t:.1f} s")
        self.truncate_after(t)
        event = apply_breakup(self.table, len(self.events), id_a, id_b, t, self.config, params)
        self.events.append(event)
        self._added.extend(event.created_ids)
        self._updated.update((event.id_a, event.id_b))
        self._added_events.append(event.event_id)
        return event

    def collide_approach(self, approach: dict, params: BreakupParams | None = None) -> BreakupEvent:
        """Force a collision of a pair from ``find_close_approaches`` at its TCA."""
        return self.trigger_breakup(approach["id_1"], approach["id_2"], approach["tca_s"], params)

    def find_close_approaches(self, t_start: float | None = None, duration_s: float = 600.0, **kwargs):
        t_start = self.clock.t if t_start is None else float(t_start)
        return screening.find_close_approaches(self.table, t_start, duration_s, **kwargs)

    def fragment_rows(self, event_id: int) -> np.ndarray:
        return np.flatnonzero(self.table.parent_event == event_id)
