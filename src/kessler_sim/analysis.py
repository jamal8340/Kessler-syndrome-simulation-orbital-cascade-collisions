"""Debris-cloud statistics at a given time (used by the demo, later by validation and Unity)."""

from __future__ import annotations

import numpy as np

from . import orbits


def cloud_stats(sim, event_id: int, t: float) -> dict:
    """How far the fragments of one breakup have spread by time ``t``."""
    tb = sim.table
    rows = sim.fragment_rows(event_id)
    rows = rows[tb.alive_mask(t)[rows]]
    if not len(rows):
        return {"event_id": event_id, "t": t, "n_alive": 0}
    r, _ = tb.states(t, rows)
    el = tb.elements.take(rows)

    radius = np.linalg.norm(r, axis=1)
    centre = r.mean(axis=0)
    centre_dir = centre / np.linalg.norm(centre) if np.linalg.norm(centre) > 0 else np.array([1.0, 0, 0])
    angle = np.degrees(np.arccos(np.clip((r / radius[:, None]) @ centre_dir, -1.0, 1.0)))

    raan_dot, _, _ = orbits.secular_rates(el)
    raan_now = np.mod(el.raan + raan_dot * (t - tb.t_ref[rows]), orbits.TWO_PI)

    per_parent = {}
    ev = sim.events[event_id]
    for pid in (ev.id_a, ev.id_b):
        sel = ev.fragment_origin_ids == pid
        ids = set(ev.fragment_ids[sel].tolist())
        mask = np.array([int(i) in ids for i in tb.ids[rows]])
        if mask.any():
            inc = np.degrees(el.i[mask])
            resultant = np.abs(np.mean(np.exp(1j * raan_now[mask])))     # circular std of RAAN
            per_parent[int(pid)] = {"n": int(mask.sum()), "i_mean_deg": float(inc.mean()),
                                    "i_std_deg": float(inc.std()),
                                    "raan_std_deg": float(np.degrees(np.sqrt(-2.0 * np.log(max(resultant, 1e-12)))))}

    return {
        "event_id": event_id, "t": float(t), "n_alive": int(len(rows)),
        "altitude_p5_km": float(np.percentile(radius, 5) - orbits.R_MEAN),
        "altitude_p95_km": float(np.percentile(radius, 95) - orbits.R_MEAN),
        "high_orbit_count": int(np.sum(orbits.apogee_radius(el) - orbits.R_MEAN > 2000.0)),
        "perigee_alt_p5_km": float(np.percentile(orbits.perigee_radius(el) - orbits.R_MEAN, 5)),
        "apogee_alt_p95_km": float(np.percentile(orbits.apogee_radius(el) - orbits.R_MEAN, 95)),
        "max_distance_from_centre_km": float(np.linalg.norm(r - centre, axis=1).max()),
        "arc_spread_p95_deg": float(np.percentile(angle, 95)),
        "per_parent": per_parent,
    }
