"""Close-approach search over a time window, using the analytic orbit model.

Same method and thresholds as ``collision_detection/close_approach.py``: sample the
window in fixed steps, find candidate pairs that could come within the threshold during
the step, and compute the time of closest approach (TCA) and miss distance assuming
linear relative motion within the step. Positions come from the orbit formula at every
step, so no integration error accumulates.

Pairs that are skipped:
* active satellites of the same constellation (coordinated formation flying);
* fragments of the same breakup (siblings), which start at the same point.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from collision_physics import fragmentation as fr

from . import orbits
from .state import NO_EVENT, ObjectTable

SCREENING_THRESHOLD_KM = 10.0
CRITICAL_CONJUNCTION_KM = 2.0
MIN_RELATIVE_SPEED_KMS = 0.5
STEP_S = 5.0
MAX_RELATIVE_SPEED_KMS = 16.0   # bounds how far two objects can close in within one step


def classify_event(miss_distance_km, hard_body_radius_km):
    if miss_distance_km <= hard_body_radius_km:
        return "PHYSICAL_COLLISION"
    if miss_distance_km <= CRITICAL_CONJUNCTION_KM:
        return "CRITICAL_CONJUNCTION"
    return "CLOSE_APPROACH"


def _candidate_pairs(r, group, radius):
    """Index pairs (i, j) closer than ``radius``, never within the same breakup group."""
    pairs = []
    groups = np.unique(group)
    trees = {}
    for g in groups:
        idx = np.flatnonzero(group == g)
        trees[g] = (idx, cKDTree(r[idx]))
    background = trees.get(NO_EVENT)
    if background is not None and len(background[0]) > 1:
        p = background[1].query_pairs(radius, output_type="ndarray")
        pairs.append(background[0][p])
    for a_pos, ga in enumerate(groups):
        for gb in groups[a_pos + 1:]:
            ia, ta = trees[ga]
            ib, tb = trees[gb]
            m = ta.sparse_distance_matrix(tb, radius, output_type="ndarray")
            if len(m):
                pairs.append(np.stack([ia[m["i"]], ib[m["j"]]], axis=1))
    if not pairs:
        return np.empty((0, 2), dtype=np.int64)
    return np.concatenate(pairs).astype(np.int64)


def find_close_approaches(table: ObjectTable, t_start: float, duration_s: float = 600.0,
                          step_s: float = STEP_S, threshold_km: float = SCREENING_THRESHOLD_KM,
                          min_rel_speed_kms: float = MIN_RELATIVE_SPEED_KMS,
                          skip_same_constellation: bool = True) -> list[dict]:
    """Closest encounter of every pair that comes within ``threshold_km`` in the window.

    Returns a list of dicts sorted by miss distance; times are absolute sim times.
    """
    radius = threshold_km + MAX_RELATIVE_SPEED_KMS * step_s
    n_steps = max(1, int(np.ceil(duration_s / step_s)))
    best: dict[tuple[int, int], dict] = {}

    for k in range(n_steps):
        t0 = t_start + k * step_s
        dt_window = min(step_s, t_start + duration_s - t0)
        rows = table.alive_rows(t0)
        if len(rows) < 2:
            continue
        r, v = table.states(t0, rows)
        # breakup group: NO_EVENT for catalogue objects, event ID for fragments/remnants
        group = table.parent_event[rows]
        pairs = _candidate_pairs(r, group, radius)
        if not len(pairs):
            continue
        i, j = pairs[:, 0], pairs[:, 1]

        if skip_same_constellation:
            ci, cj = table.constellation[rows[i]], table.constellation[rows[j]]
            keep = ~((ci == cj) & (ci != ""))
            i, j = i[keep], j[keep]

        dr = r[i] - r[j]
        dv = v[i] - v[j]
        dv2 = np.einsum("ij,ij->i", dv, dv)
        t_ca = np.where(dv2 > 1e-12, -np.einsum("ij,ij->i", dr, dv) / np.maximum(dv2, 1e-12), 0.0)
        t_ca = np.clip(t_ca, 0.0, dt_window)
        d_min = np.linalg.norm(dr + dv * t_ca[:, None], axis=1)
        rel = np.sqrt(dv2)
        hit = (d_min <= threshold_km) & (rel >= min_rel_speed_kms)

        for a, b, tc, d, vr in zip(rows[i[hit]], rows[j[hit]], t_ca[hit], d_min[hit], rel[hit]):
            id_a, id_b = int(table.ids[a]), int(table.ids[b])
            if id_a > id_b:
                id_a, id_b, a, b = id_b, id_a, b, a
            key = (id_a, id_b)
            if key in best and best[key]["miss_distance_km"] <= d:
                continue
            best[key] = {"_rows": (a, b), "tca_s": t0 + float(tc),
                         "miss_distance_km": float(d), "relative_speed_kms": float(vr)}

    results = []
    for (id_a, id_b), rec in best.items():
        a, b = rec.pop("_rows")
        (ra, rb), _ = table.states(rec["tca_s"], [a, b])
        midpoint = 0.5 * (ra + rb)
        hbr_km = (table.radius_m[a] + table.radius_m[b]) / 1000.0
        _, emr, cata = fr.collision_energetics(table.mass_kg[a], table.mass_kg[b], rec["relative_speed_kms"])
        results.append({
            "id_1": id_a, "name_1": str(table.names[a]), "category_1": str(table.category[a]),
            "id_2": id_b, "name_2": str(table.names[b]), "category_2": str(table.category[b]),
            "event_type": classify_event(rec["miss_distance_km"], hbr_km),
            **rec,
            "altitude_km": float(np.linalg.norm(midpoint) - orbits.R_MEAN),
            "position_km": midpoint.tolist(),
            "hard_body_radius_km": float(hbr_km),
            "emr_jg": float(emr), "would_be_catastrophic": bool(cata),
            "debris_involved": str(table.category[a]) != "active" or str(table.category[b]) != "active",
        })
    results.sort(key=lambda e: e["miss_distance_km"])
    return results
