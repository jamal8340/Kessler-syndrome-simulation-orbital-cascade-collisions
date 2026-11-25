"""Turning a collision into new orbiting objects.

Uses the NASA Standard Breakup Model from ``collision_physics.fragmentation`` for the
fragment count, sizes, masses and velocity kicks (dv), then:

* Catastrophic collision: each fragment is assigned to one of the two parents (largest
  fragments first, always to the parent with the most unassigned mass), and its velocity
  is the parent's velocity plus its dv. The dvs are re-centred per parent, so each
  parent's momentum is conserved exactly and its debris keeps following the parent's
  orbit (two clouds, as observed for real collisions).
* Non-catastrophic collision: fragments are ejected around the centre-of-mass velocity
  (practically the target's velocity) and the target survives as a new "remnant" object
  carrying the remaining mass.
* All fragments start at the heavier object's position. For a pair taken from the
  close-approach list at its TCA, the two objects are within a few km of each other;
  for an arbitrary "what if" pair, ``separation_km`` in the event shows how artificial
  the setup is.
* Fragments whose perigee is below the re-entry altitude, or that are not on a bound
  orbit, are counted but not added (crude drag model).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from collision_physics import fragmentation as fr

from . import orbits
from .state import FRAGMENT_CATEGORY, ObjectTable

FRAGMENT_ID_BASE = 1_000_000
FRAGMENT_IDS_PER_EVENT = 100_000


@dataclass
class BreakupConfig:
    seed: int = 42
    l_min_m: float = fr.L_MIN_DEFAULT_M
    max_fragments: int = 5000
    reentry_altitude_km: float = 150.0


@dataclass
class BreakupParams:
    """Optional per-breakup overrides ("what if" settings)."""

    mass_a_kg: float | None = None
    mass_b_kg: float | None = None
    impact_speed_kms: float | None = None
    l_min_m: float | None = None
    max_fragments: int | None = None
    seed: int | None = None


@dataclass
class BreakupEvent:
    event_id: int
    t: float
    id_a: int
    id_b: int
    name_a: str
    name_b: str
    mass_a_kg: float
    mass_b_kg: float
    position_km: np.ndarray
    separation_km: float
    v_rel_kms: float
    v_cm_kms: np.ndarray
    catastrophic: bool
    energy_j: float
    emr_jg: float
    l_min_m: float
    n_expected: int
    seed: int
    fragment_ids: np.ndarray = field(repr=False)
    fragment_origin_ids: np.ndarray = field(repr=False)   # parent each fragment belongs to
    fragment_dv_kms: np.ndarray = field(repr=False)       # velocity relative to its parent / CM
    fragment_mass_kg: np.ndarray = field(repr=False)
    fragment_lc_m: np.ndarray = field(repr=False)
    survivor_id: int | None
    tracked_mass_kg: float
    untracked_mass_kg: float
    reentered_count: int
    escaped_count: int
    lost_mass_kg: float          # mass of fragments that re-entered or escaped

    @property
    def n_fragments(self) -> int:
        return len(self.fragment_ids)

    @property
    def created_ids(self) -> list[int]:
        ids = self.fragment_ids.tolist()
        return ids + ([self.survivor_id] if self.survivor_id is not None else [])

    def summary(self) -> dict:
        """JSON-friendly summary (no per-fragment arrays)."""
        largest = int(np.argmax(self.fragment_mass_kg)) if self.n_fragments else None
        per_parent = {
            int(pid): int(np.sum(self.fragment_origin_ids == pid)) for pid in (self.id_a, self.id_b)
        }
        return {
            "event_id": self.event_id, "t": self.t,
            "id_a": self.id_a, "id_b": self.id_b, "name_a": self.name_a, "name_b": self.name_b,
            "mass_a_kg": self.mass_a_kg, "mass_b_kg": self.mass_b_kg,
            "position_km": self.position_km.tolist(),
            "altitude_km": float(np.linalg.norm(self.position_km) - orbits.R_MEAN),
            "separation_km": self.separation_km, "v_rel_kms": self.v_rel_kms,
            "catastrophic": self.catastrophic, "energy_j": self.energy_j, "emr_jg": self.emr_jg,
            "l_min_m": self.l_min_m, "n_expected": self.n_expected, "n_fragments": self.n_fragments,
            "fragments_per_parent": per_parent,
            "survivor_id": self.survivor_id,
            "tracked_mass_kg": self.tracked_mass_kg, "untracked_mass_kg": self.untracked_mass_kg,
            "reentered_count": self.reentered_count, "escaped_count": self.escaped_count,
            "lost_mass_kg": self.lost_mass_kg,
            "largest_fragment_mass_kg": float(self.fragment_mass_kg[largest]) if largest is not None else None,
            "largest_fragment_lc_m": float(self.fragment_lc_m[largest]) if largest is not None else None,
            "seed": self.seed,
        }


def fragment_ids_for(event_id: int, count: int) -> np.ndarray:
    if count >= FRAGMENT_IDS_PER_EVENT:
        raise ValueError(f"too many fragments for one event ({count})")
    return FRAGMENT_ID_BASE + event_id * FRAGMENT_IDS_PER_EVENT + np.arange(count, dtype=np.int64)


def survivor_id_for(event_id: int) -> int:
    return FRAGMENT_ID_BASE + event_id * FRAGMENT_IDS_PER_EVENT + FRAGMENT_IDS_PER_EVENT - 1


def capped_l_min(m_a, m_b, v_rel_kms, l_min_m, max_fragments):
    """Raise the minimum fragment size so the expected count does not exceed the cap."""
    m_t, m_p = max(m_a, m_b), min(m_a, m_b)
    _, _, cata = fr.collision_energetics(m_t, m_p, v_rel_kms)
    m_inv = fr.mass_involved(m_t, m_p, v_rel_kms, cata)
    if m_inv <= 0 or fr.expected_fragment_count(m_inv, l_min_m) <= max_fragments:
        return l_min_m
    return (max_fragments / (fr.N_COEFF * m_inv ** fr.N_MASS_EXP)) ** (-1.0 / fr.N_SIZE_EXP)


def _assign_to_parents(masses, budget_a, budget_b):
    """Greedy split: largest fragment first, to the parent with more unassigned mass."""
    remaining = np.array([budget_a, budget_b], dtype=float)
    owner = np.empty(len(masses), dtype=np.int64)
    for k in np.argsort(masses)[::-1]:
        side = int(remaining[1] > remaining[0])
        owner[k] = side
        remaining[side] -= masses[k]
    return owner, remaining


def apply_breakup(table: ObjectTable, event_id: int, id_a: int, id_b: int, t: float,
                  config: BreakupConfig, params: BreakupParams | None = None) -> BreakupEvent:
    params = params or BreakupParams()
    if int(id_a) == int(id_b):
        raise ValueError("an object cannot collide with itself")
    row_a, row_b = table.row(id_a), table.row(id_b)
    for oid, row in ((id_a, row_a), (id_b, row_b)):
        if not table.is_alive(row, t):
            raise ValueError(f"object {oid} does not exist at t = {t:.1f} s")

    (r_a, r_b), (v_a, v_b) = table.states(t, [row_a, row_b])
    m_a = float(params.mass_a_kg if params.mass_a_kg is not None else table.mass_kg[row_a])
    m_b = float(params.mass_b_kg if params.mass_b_kg is not None else table.mass_kg[row_b])
    if m_a <= 0 or m_b <= 0:
        raise ValueError("masses must be positive")

    # heavier object = target; its position is the collision point
    a_is_target = m_a >= m_b
    r_t = r_a if a_is_target else r_b
    v_t, v_p = (v_a, v_b) if a_is_target else (v_b, v_a)
    if params.impact_speed_kms is not None:
        rel = v_p - v_t
        norm = np.linalg.norm(rel)
        direction = rel / norm if norm > 0 else np.array([1.0, 0.0, 0.0])
        v_p = v_t + direction * float(params.impact_speed_kms)
        if a_is_target:
            v_b = v_p
        else:
            v_a = v_p
    v_rel = float(np.linalg.norm(v_a - v_b))

    l_min = float(params.l_min_m if params.l_min_m is not None else config.l_min_m)
    cap = int(params.max_fragments if params.max_fragments is not None else config.max_fragments)
    l_min = capped_l_min(m_a, m_b, v_rel, l_min, min(cap, FRAGMENT_IDS_PER_EVENT - 1))

    seed = int(params.seed) if params.seed is not None else int(
        np.random.SeedSequence([config.seed, event_id]).generate_state(1)[0])
    body_a = fr.Body(m_a, float(table.radius_m[row_a]), tuple(v_a))
    body_b = fr.Body(m_b, float(table.radius_m[row_b]), tuple(v_b))
    res = fr.simulate_collision(body_a, body_b, l_min_m=l_min, seed=seed)

    n = res.n_fragments
    v_cm = res.v_cm_kms
    dv = res.velocity_kms - v_cm
    origin = np.full(n, int(id_a if a_is_target else id_b), dtype=np.int64)
    untracked = res.remnant_mass_kg
    survivor = None

    if res.catastrophic and n:
        owner, remaining = _assign_to_parents(res.mass_kg, m_a, m_b)
        velocity = np.empty_like(dv)
        for side, (oid, v_parent) in enumerate(((id_a, v_a), (id_b, v_b))):
            sel = owner == side
            if sel.any():
                dv[sel] -= (res.mass_kg[sel, None] * dv[sel]).sum(axis=0) / res.mass_kg[sel].sum()
                velocity[sel] = v_parent + dv[sel]
                origin[sel] = oid
        untracked = float(np.clip(remaining, 0.0, None).sum())
    else:
        velocity = res.velocity_kms

    # fragments -> orbits; drop the ones that re-enter or escape
    positions = np.broadcast_to(r_t, (n, 3))
    el, bound = orbits.state_to_elements(positions, velocity) if n else (orbits.Elements.empty(), np.empty(0, bool))
    perigee_ok = orbits.perigee_radius(el) > orbits.R_MEAN + config.reentry_altitude_km if n else bound
    keep = bound & perigee_ok
    escaped = int(np.sum(~bound))
    reentered = int(np.sum(bound & ~perigee_ok))

    kept = np.flatnonzero(keep)
    frag_ids = fragment_ids_for(event_id, len(kept))
    if len(kept):
        name_of = {int(id_a): str(table.names[row_a]), int(id_b): str(table.names[row_b])}
        table.add(
            frag_ids,
            names=[f"{name_of[int(o)]} DEB #{k + 1}" for k, o in enumerate(origin[kept])],
            category=FRAGMENT_CATEGORY,
            mass_kg=res.mass_kg[kept], radius_m=res.lc_m[kept] / 2.0,
            elements=el.take(kept), t_ref=t, t_created=t,
            parent_event=event_id, parent_ids=(int(id_a), int(id_b)),
            am_m2_per_kg=res.am_m2_per_kg[kept],
        )

    if not res.catastrophic:
        target_row = row_a if a_is_target else row_b
        survivor = survivor_id_for(event_id)
        s_el, _ = orbits.state_to_elements(r_t, v_cm)
        table.add(
            [survivor], names=f"{table.names[target_row]} [REMNANT]",
            category=table.category[target_row], mass_kg=res.remnant_mass_kg,
            radius_m=table.radius_m[target_row], elements=s_el, t_ref=t, t_created=t,
            parent_event=event_id, parent_ids=(int(id_a), int(id_b)),
        )
        untracked = 0.0

    table.destroy([row_a, row_b], t)

    return BreakupEvent(
        event_id=event_id, t=float(t), id_a=int(id_a), id_b=int(id_b),
        name_a=str(table.names[row_a]), name_b=str(table.names[row_b]),
        mass_a_kg=m_a, mass_b_kg=m_b,
        position_km=np.array(r_t), separation_km=float(np.linalg.norm(r_a - r_b)),
        v_rel_kms=v_rel, v_cm_kms=np.array(v_cm),
        catastrophic=bool(res.catastrophic), energy_j=float(res.energy_j), emr_jg=float(res.emr_jg),
        l_min_m=l_min, n_expected=int(res.n_expected), seed=seed,
        fragment_ids=frag_ids, fragment_origin_ids=origin[kept],
        fragment_dv_kms=dv[kept], fragment_mass_kg=res.mass_kg[kept], fragment_lc_m=res.lc_m[kept],
        survivor_id=survivor,
        tracked_mass_kg=float(res.mass_kg[kept].sum()),
        untracked_mass_kg=float(untracked),
        reentered_count=reentered, escaped_count=escaped,
        lost_mass_kg=float(res.mass_kg[~keep].sum()) if n else 0.0,
    )
