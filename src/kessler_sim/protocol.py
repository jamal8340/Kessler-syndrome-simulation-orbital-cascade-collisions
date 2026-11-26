"""Messages between the backend and Unity (see TRANSPORT_ARCHITECTURE.md).

Pure functions: they turn engine data into messages and back, with no networking,
so they can be tested on their own.

* JSON messages never contain infinity or NaN: such values become ``null``.
* State messages are binary (little-endian): float64 t, float64 earth_rotation_rad,
  uint32 catalogue_version, uint32 count, then count x 3 float32 positions [km]
  in catalogue order, NaN for objects that do not exist at t.
"""

from __future__ import annotations

import json
import math
import struct

import numpy as np

from .earth import earth_rotation_at
from .engine import CatalogueChanges, Simulation

STATE_HEADER = struct.Struct("<ddII")
FRAME = "TEME"
UNITS = {"time": "s since epoch", "position": "km", "velocity": "km/s", "mass": "kg",
         "radius": "m", "angle": "rad"}


# -- JSON -------------------------------------------------------------------------

def _clean(value):
    """Make a value JSON-safe: numpy -> Python, inf/NaN -> None, tuples -> lists."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return _clean(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        return value if math.isfinite(value) else None
    return value


def dumps(message: dict) -> str:
    """JSON text of a message; raises if anything non-finite slipped through."""
    return json.dumps(_clean(message), allow_nan=False, separators=(",", ":"))


def message(kind: str, request_id=None, **fields) -> dict:
    msg = {"type": kind}
    if request_id is not None:
        msg["request_id"] = request_id
    msg.update(fields)
    return msg


# -- catalogue --------------------------------------------------------------------

def catalogue_entry(sim: Simulation, row: int) -> dict:
    tb = sim.table
    return _clean({
        "id": int(tb.ids[row]),
        "name": str(tb.names[row]),
        "category": str(tb.category[row]),
        "mass_kg": float(tb.mass_kg[row]),
        "radius_m": float(tb.radius_m[row]),
        "t_created": float(tb.t_created[row]),        # -inf -> null (exists from the start)
        "t_destroyed": float(tb.t_destroyed[row]),    # +inf -> null (not destroyed)
        "parent_event": None if tb.parent_event[row] < 0 else int(tb.parent_event[row]),
        "parent_ids": list(tb.parent_ids[row]),
    })


def catalogue(sim: Simulation) -> list[dict]:
    return [catalogue_entry(sim, row) for row in range(len(sim.table))]


def catalogue_update(sim: Simulation, changes: CatalogueChanges) -> dict:
    tb = sim.table
    return message(
        "catalogue_update",
        catalogue_version=tb.version,
        removed=changes.removed,
        updated=[catalogue_entry(sim, tb.row(i)) for i in changes.updated],
        added=[catalogue_entry(sim, tb.row(i)) for i in changes.added],
        count=len(tb),
    )


# -- events -------------------------------------------------------------------------

def breakup_event(event, request_id=None) -> dict:
    return message(
        "breakup_event", request_id,
        **event.summary(),
        remnant_id=event.survivor_id,
        v_cm_kms=event.v_cm_kms,
        fragment_ids=event.fragment_ids,
        fragment_parent_ids=event.fragment_origin_ids,
        fragment_dv_kms=event.fragment_dv_kms,
        fragment_mass_kg=event.fragment_mass_kg,
        fragment_lc_m=event.fragment_lc_m,
    )


def events_changed(removed_event_ids, request_id=None) -> dict:
    return message("events_changed", request_id, removed=list(removed_event_ids))


def clock_status(sim: Simulation, request_id=None) -> dict:
    return message("clock_status", request_id, **sim.clock.status(),
                   utc=sim.utc_at(sim.clock.t).isoformat())


def session(sim: Simulation, kind: str = "hello", request_id=None) -> dict:
    """Everything a client needs to (re)start: sent as the hello reply and after load_scenario."""
    return message(
        kind, request_id,
        scenario={"name": sim.name, "description": getattr(sim, "description", "")},
        epoch_utc=sim.epoch_utc.isoformat(),
        epoch_jd=sim.epoch_jd,
        frame=FRAME,
        units=UNITS,
        clock=sim.clock.status(),
        catalogue_version=sim.table.version,
        catalogue=catalogue(sim),
        events=[breakup_event(ev) for ev in sim.events],
    )


def error(reason: str, request_id=None) -> dict:
    return message("error", request_id, reason=reason)


# -- binary state -------------------------------------------------------------------

def state_positions(sim: Simulation, t: float) -> np.ndarray:
    """(N, 3) float32 positions of every catalogue entry at t, NaN where it does not exist."""
    tb = sim.table
    rows = tb.alive_rows(t)
    r = np.full((len(tb), 3), np.nan, dtype=np.float32)
    r[rows], _ = tb.states(t, rows)
    return r


def pack_state(t: float, earth_rotation_rad: float, version: int, positions: np.ndarray) -> bytes:
    positions = np.ascontiguousarray(positions, dtype="<f4")
    return STATE_HEADER.pack(float(t), float(earth_rotation_rad), int(version), len(positions)) + positions.tobytes()


def state_message(sim: Simulation, t: float | None = None) -> bytes:
    t = sim.clock.t if t is None else float(t)
    return pack_state(t, earth_rotation_at(sim, t), sim.table.version, state_positions(sim, t))


def unpack_state(data: bytes) -> dict:
    t, angle, version, count = STATE_HEADER.unpack_from(data)
    positions = np.frombuffer(data, dtype="<f4", count=3 * count, offset=STATE_HEADER.size).reshape(count, 3)
    return {"t": t, "earth_rotation_rad": angle, "catalogue_version": version, "positions": positions}
