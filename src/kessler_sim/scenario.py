"""Scenario files: what to load and which breakups to schedule.

A scenario is a JSON file (see ``scenarios/*.json``)::

    {
      "name": "...", "description": "...",
      "catalogue": true,                     # load the shared dataset (its epoch is t = 0)
      "epoch_utc": "2009-02-10T16:46:00Z",   # required when "catalogue" is false
      "seed": 42, "l_min_m": 0.1, "max_fragments": 5000, "reentry_altitude_km": 150,
      "objects": [                           # extra objects, defined by mean elements
        {"id": 24946, "name": "IRIDIUM 33", "category": "active", "mass_kg": 560,
         "radius_m": 1.2, "t_ref_s": 600,
         "elements": {"a_km": 7160, "e": 0, "i_deg": 86.4, "raan_deg": 348.49,
                      "argp_deg": 0, "mean_anomaly_deg": 72.86}}
      ],
      "breakups": [
        {"t_s": 600, "id_a": 24946, "id_b": 22675, "params": {"mass_a_kg": 560}},
        {"closest_approach": {"t_start_s": 0, "duration_s": 600}}
      ]
    }

Instead of ``elements``, an object may be given by ``"state": {"r_km": [x, y, z], "v_kms": [vx, vy, vz]}``
at ``t_ref_s`` (any bound orbit, also a suborbital one), and by ``"t_created_s"`` it only exists
from that time on (e.g. an interceptor's launch).

``closest_approach`` searches the window and collides the closest pair at its TCA.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np

from . import orbits
from .breakup import BreakupConfig, BreakupParams
from .catalogue import datetime_to_jd, load_catalogue
from .engine import Simulation
from .state import ObjectTable

SCENARIO_DIR = Path(__file__).resolve().parent / "scenarios"


def available_scenarios() -> list[str]:
    return sorted(p.stem for p in SCENARIO_DIR.glob("*.json"))


def resolve(name_or_path) -> Path:
    path = Path(name_or_path)
    if path.suffix == ".json" and path.is_file():
        return path
    candidate = SCENARIO_DIR / f"{name_or_path}.json"
    if candidate.is_file():
        return candidate
    raise FileNotFoundError(f"scenario '{name_or_path}' not found; available: {available_scenarios()}")


def _elements_from_spec(spec: dict) -> orbits.Elements:
    if "state" in spec:
        # position [km] and velocity [km/s] at t_ref_s (e.g. a suborbital interceptor)
        st = spec["state"]
        el, bound = orbits.state_to_elements(np.array([st["r_km"]], float), np.array([st["v_kms"]], float))
        if not bound[0]:
            raise ValueError(f"object {spec.get('id')}: state is not on a bound orbit")
        return el
    el = spec["elements"]
    return orbits.Elements(
        np.array([el["a_km"]], float), np.array([el.get("e", 0.0)], float),
        np.radians([el["i_deg"]]), np.radians([el.get("raan_deg", 0.0)]),
        np.radians([el.get("argp_deg", 0.0)]), np.radians([el.get("mean_anomaly_deg", 0.0)]),
    )


def load_scenario(name_or_path) -> Simulation:
    path = resolve(name_or_path)
    spec = json.loads(path.read_text(encoding="utf-8"))

    config = BreakupConfig(
        seed=int(spec.get("seed", BreakupConfig.seed)),
        l_min_m=float(spec.get("l_min_m", BreakupConfig.l_min_m)),
        max_fragments=int(spec.get("max_fragments", BreakupConfig.max_fragments)),
        reentry_altitude_km=float(spec.get("reentry_altitude_km", BreakupConfig.reentry_altitude_km)),
    )

    table = ObjectTable()
    report = None
    if spec.get("catalogue", False):
        _, epoch_jd, report = load_catalogue(table=table)
    else:
        epoch = datetime.fromisoformat(spec["epoch_utc"].replace("Z", "+00:00"))
        epoch_jd = datetime_to_jd(epoch)

    for obj in spec.get("objects", []):
        t_created = obj.get("t_created_s")          # e.g. launch time; default: exists from the start
        table.add([obj["id"]], names=obj["name"], category=obj.get("category", "active"),
                  mass_kg=obj["mass_kg"], radius_m=obj["radius_m"],
                  elements=_elements_from_spec(obj), t_ref=float(obj.get("t_ref_s", 0.0)),
                  t_created=-np.inf if t_created is None else float(t_created))

    sim = Simulation(table, epoch_jd, config, name=spec.get("name", path.stem))
    sim.description = spec.get("description", "")
    sim.load_report = report

    for b in spec.get("breakups", []):
        params = BreakupParams(**b.get("params", {}))
        if "closest_approach" in b:
            window = b["closest_approach"]
            approaches = sim.find_close_approaches(window.get("t_start_s", 0.0), window.get("duration_s", 600.0))
            if not approaches:
                raise ValueError(f"scenario {path.stem}: no close approach found in {window}")
            sim.collide_approach(approaches[0], params)
        else:
            sim.trigger_breakup(b["id_a"], b["id_b"], float(b["t_s"]), params)
    return sim
