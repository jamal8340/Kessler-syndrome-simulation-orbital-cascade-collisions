"""Simulation engine: the single entry point the server (and scripts, tests) use.

It owns the object table, the clock and the event log. Positions are computed on
demand for any time t, so play/pause/seek are trivial and runs are reproducible:
the state is fully determined by the dataset, the seed and the list of breakups.

Rewriting history: triggering a breakup at time t discards every breakup that
happened after t (their fragments are removed and their parents come back).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import orbits, screening
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


class Simulation:
    def __init__(self, table: ObjectTable | None = None, epoch_jd: float = 2451545.0,
                 config: BreakupConfig | None = None, name: str = ""):
        self.table = table if table is not None else ObjectTable()
        self.epoch_jd = float(epoch_jd)
        self.config = config or BreakupConfig()
        self.name = name
        self.clock = Clock()
        self.events: list[BreakupEvent] = []

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
        el = tb.elements.take([row])
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

    # -- events ------------------------------------------------------------------

    def truncate_after(self, t: float) -> list[BreakupEvent]:
        """Remove every breakup after ``t`` (newest first) and restore its parents."""
        removed = []
        while self.events and self.events[-1].t > t:
            ev = self.events.pop()
            created = [self.table.row(i) for i in ev.created_ids if self.table.has(i)]
            parents = [self.table.row(ev.id_a), self.table.row(ev.id_b)]
            self.table.revive(parents)
            if created:
                self.table.remove_rows(created)
            removed.append(ev)
        return removed

    def trigger_breakup(self, id_a: int, id_b: int, t: float | None = None,
                        params: BreakupParams | None = None) -> BreakupEvent:
        t = self.clock.t if t is None else float(t)
        # validate before rewriting history, so an invalid request never discards events
        if int(id_a) == int(id_b):
            raise ValueError("an object cannot collide with itself")
        for oid in (id_a, id_b):
            if not self.table.is_alive(self.table.row(oid), t):
                raise ValueError(f"object {oid} does not exist at t = {t:.1f} s")
        self.truncate_after(t)
        event = apply_breakup(self.table, len(self.events), id_a, id_b, t, self.config, params)
        self.events.append(event)
        return event

    def collide_approach(self, approach: dict, params: BreakupParams | None = None) -> BreakupEvent:
        """Force a collision of a pair from ``find_close_approaches`` at its TCA."""
        return self.trigger_breakup(approach["id_1"], approach["id_2"], approach["tca_s"], params)

    def find_close_approaches(self, t_start: float | None = None, duration_s: float = 600.0, **kwargs):
        t_start = self.clock.t if t_start is None else float(t_start)
        return screening.find_close_approaches(self.table, t_start, duration_s, **kwargs)

    def fragment_rows(self, event_id: int) -> np.ndarray:
        return np.flatnonzero(self.table.parent_event == event_id)
