"""Object table: every object in the simulation (catalogue objects and fragments).

Each object has an orbit (mean elements referenced to ``t_ref``) and a lifespan
[t_created, t_destroyed). The state of the simulation at time t is simply "all objects
alive at t, at their positions at t", so seeking backwards before a breakup shows the
parents again.
"""

from __future__ import annotations

import numpy as np

from . import drag, orbits
from .orbit_history import OrbitHistory

CONSTELLATION_PREFIXES = ("STARLINK", "ONEWEB", "IRIDIUM", "FLOCK")

FRAGMENT_CATEGORY = "fragment"
NO_EVENT = -1


def constellation_of(name: str, category: str) -> str:
    """Constellation tag used to skip formation-flying pairs ('' = none)."""
    if category != "active":
        return ""
    upper = name.upper()
    for prefix in CONSTELLATION_PREFIXES:
        if upper.startswith(prefix):
            return prefix
    return ""


class ObjectTable:
    """Struct-of-arrays table of objects. Rows are positional; IDs are stable."""

    def __init__(self):
        self.ids = np.empty(0, dtype=np.int64)
        self.names = np.empty(0, dtype=object)
        self.category = np.empty(0, dtype=object)
        self.constellation = np.empty(0, dtype=object)
        self.mass_kg = np.empty(0)
        self.radius_m = np.empty(0)
        self.am_m2_per_kg = np.empty(0)          # NaN for catalogue objects
        self.elements = orbits.Elements.empty()
        self.t_ref = np.empty(0)
        self.t_created = np.empty(0)
        self.t_destroyed = np.empty(0)
        self.parent_event = np.empty(0, dtype=np.int64)
        self.parent_ids = np.empty(0, dtype=object)  # tuple of parent IDs ( () for catalogue )
        self._row_of: dict[int, int] = {}
        self.orbit_history = OrbitHistory()
        self.version = 0                          # bumped on every structural change

    _COLUMNS = ("ids", "names", "category", "constellation", "mass_kg", "radius_m", "am_m2_per_kg",
                "t_ref", "t_created", "t_destroyed", "parent_event", "parent_ids")

    def __len__(self):
        return len(self.ids)

    # -- building ------------------------------------------------------------

    def add(self, ids, names, category, mass_kg, radius_m, elements: orbits.Elements, t_ref,
            t_created=-np.inf, parent_event=NO_EVENT, parent_ids=(), am_m2_per_kg=np.nan):
        """Append objects; scalar arguments are broadcast. Returns the new rows."""
        ids = np.asarray(ids, dtype=np.int64)
        n = len(ids)
        if len(elements) != n:
            raise ValueError("elements length does not match ids")
        duplicates = [int(i) for i in ids if int(i) in self._row_of]
        if duplicates or len(set(ids.tolist())) != n:
            raise ValueError(f"duplicate object IDs: {duplicates[:5]}")

        def col(x, dtype=float):
            if dtype is object:
                out = np.empty(n, dtype=object)
                if isinstance(x, (str, tuple)) or x is None:
                    for k in range(n):
                        out[k] = x
                else:
                    values = list(x)
                    for k in range(n):
                        out[k] = values[k]
                return out
            return np.broadcast_to(np.asarray(x, dtype=dtype), (n,)).copy()

        names = col(names, object)
        category = col(category, object)
        start = len(self.ids)
        self.ids = np.concatenate([self.ids, ids])
        self.names = np.concatenate([self.names, names])
        self.category = np.concatenate([self.category, category])
        self.constellation = np.concatenate(
            [self.constellation, col([constellation_of(nm, c) for nm, c in zip(names, category)], object)])
        self.mass_kg = np.concatenate([self.mass_kg, col(mass_kg)])
        self.radius_m = np.concatenate([self.radius_m, col(radius_m)])
        self.am_m2_per_kg = np.concatenate([self.am_m2_per_kg, col(am_m2_per_kg)])
        self.elements = orbits.Elements.concat([self.elements, elements])
        self.t_ref = np.concatenate([self.t_ref, col(t_ref)])
        self.t_created = np.concatenate([self.t_created, col(t_created)])
        self.t_destroyed = np.concatenate([self.t_destroyed, np.full(n, np.inf)])
        self.parent_event = np.concatenate([self.parent_event, col(parent_event, np.int64)])
        self.parent_ids = np.concatenate([self.parent_ids, col(tuple(parent_ids), object)])
        rows = np.arange(start, start + n)
        self._row_of.update(zip(ids.tolist(), rows.tolist()))
        self.orbit_history.add(n)
        self.version += 1
        return rows

    def remove_rows(self, rows):
        """Delete rows entirely (used when history is rewritten)."""
        keep = np.ones(len(self.ids), dtype=bool)
        keep[np.asarray(rows, dtype=np.int64)] = False
        self.orbit_history.remove(self.ids, keep)
        for name in self._COLUMNS:
            setattr(self, name, getattr(self, name)[keep])
        self.elements = self.elements.take(keep)
        self._row_of = {int(i): r for r, i in enumerate(self.ids)}
        self.version += 1

    def destroy(self, rows, t):
        self.t_destroyed[np.asarray(rows)] = t
        self.version += 1

    def revive(self, rows):
        self.t_destroyed[np.asarray(rows)] = np.inf
        self.version += 1

    # -- queries -------------------------------------------------------------

    def row(self, object_id: int) -> int:
        try:
            return self._row_of[int(object_id)]
        except KeyError:
            raise KeyError(f"unknown object ID {object_id}") from None

    def has(self, object_id: int) -> bool:
        return int(object_id) in self._row_of

    def alive_mask(self, t: float) -> np.ndarray:
        return (self.t_created <= t) & (t < self.t_destroyed)

    def alive_rows(self, t: float) -> np.ndarray:
        return np.flatnonzero(self.alive_mask(t))

    def is_alive(self, row: int, t: float) -> bool:
        return bool(self.t_created[row] <= t < self.t_destroyed[row])

    def elements_at(self, t: float, rows=None):
        """Mean elements at t, including any valid drag segment (angles at t)."""
        rows = np.arange(len(self.ids)) if rows is None else np.asarray(rows, dtype=np.int64)
        el = orbits.rebase(self.elements.take(rows), t - self.t_ref[rows])
        segments = self.orbit_history.at(self.ids, rows, t)
        use = np.isfinite(segments[:, drag.START])
        if use.any():
            current = drag.segment_elements(segments[use], t)
            for field in orbits.Elements.FIELDS:
                getattr(el, field)[use] = getattr(current, field)
        return el

    def states(self, t: float, rows=None):
        """Position and velocity (n, 3) [km, km/s] of the given rows at time t."""
        if rows is None:
            rows = np.arange(len(self.ids))
        rows = np.asarray(rows, dtype=np.int64)
        if not np.isfinite(self.orbit_history.latest[rows, drag.START]).any():
            el = self.elements.take(rows)
            return orbits.elements_to_state(el, t - self.t_ref[rows])
        return orbits.elements_to_state(self.elements_at(t, rows), 0.)
