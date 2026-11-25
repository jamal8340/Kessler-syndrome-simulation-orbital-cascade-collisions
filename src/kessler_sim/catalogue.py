"""Loading the shared catalogue dataset into an ObjectTable.

The object list (IDs, names, categories, estimated masses and sizes) comes from the
shared dataset ``data/processed/satellites_state.json``; the epoch from its summary.
Orbits come from the raw TLE files: SGP4 is evaluated at the epoch and its *mean*
elements at that moment are taken (they already include drag up to the epoch). SGP4 is
evaluated once more one day later to estimate each object's mean-motion decay, which is
the only drag effect in the model.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from sgp4.api import Satrec

from . import orbits
from .state import ObjectTable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATASET = PROJECT_ROOT / "data" / "processed" / "satellites_state.json"
DEFAULT_SUMMARY = PROJECT_ROOT / "data" / "metadata" / "simulation_summary.json"
DEFAULT_TLE_DIR = PROJECT_ROOT / "data" / "raw" / "tle"

# same order as data_pipeline/process_TLE_to_json.py: the first occurrence of an ID wins
TLE_FILES_ORDER = ("cosmos-2251-debris.tle", "iridium-33-debris.tle", "fengyun-1c-debris.tle",
                   "cosmos-1408-debris.tle", "active.tle")

DECAY_FIT_DAYS = 1.0
_JD_UNIX_EPOCH = 2440587.5


def jd_to_datetime(jd: float) -> datetime:
    # a float64 JD resolves ~40 us; round to whole milliseconds so 16:46:00 stays 16:46:00
    ms = round((jd - _JD_UNIX_EPOCH) * 86400.0 * 1000.0)
    return datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(milliseconds=ms)


def datetime_to_jd(dt: datetime) -> float:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return _JD_UNIX_EPOCH + (dt - datetime(1970, 1, 1, tzinfo=timezone.utc)).total_seconds() / 86400.0


def _split_jd(jd: float):
    whole = np.floor(jd - 0.5) + 0.5
    return float(whole), float(jd - whole)


@dataclass
class CatalogueLoadReport:
    dataset_objects: int = 0
    loaded: int = 0
    missing_tle: list[int] = field(default_factory=list)
    sgp4_failed: list[int] = field(default_factory=list)


def read_tle_records(tle_dir: Path) -> dict[int, Satrec]:
    records: dict[int, Satrec] = {}
    for fname in TLE_FILES_ORDER:
        path = tle_dir / fname
        if not path.is_file():
            continue
        lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
        for k in range(len(lines) - 1):
            if lines[k].startswith("1 ") and lines[k + 1].startswith("2 "):
                sat = Satrec.twoline2rv(lines[k], lines[k + 1])
                records.setdefault(int(sat.satnum), sat)
    return records


def mean_elements_at(sat: Satrec, jd: float):
    """SGP4 mean elements (a, e, i, raan, argp, M) and mean motion [rad/s] at ``jd``."""
    err, _, _ = sat.sgp4(*_split_jd(jd))
    if err != 0:
        return None
    n = sat.nm / 60.0
    a = (orbits.MU / n ** 2) ** (1.0 / 3.0)
    return (a, sat.em, sat.im, sat.Om, sat.om, sat.mm), n


def load_catalogue(dataset: Path = DEFAULT_DATASET, summary: Path = DEFAULT_SUMMARY,
                   tle_dir: Path = DEFAULT_TLE_DIR, table: ObjectTable | None = None,
                   limit: int | None = None):
    """Returns (table, epoch_jd, report). Sim time t = 0 is the dataset epoch."""
    objects = json.loads(Path(dataset).read_text(encoding="utf-8"))
    if limit:
        objects = objects[:limit]
    epoch_jd = float(json.loads(Path(summary).read_text(encoding="utf-8"))["propagation_epoch_jd"])
    tles = read_tle_records(Path(tle_dir))

    report = CatalogueLoadReport(dataset_objects=len(objects))
    rows, els = [], []
    for obj in objects:
        norad = int(obj["norad_id"])
        if norad not in tles:
            report.missing_tle.append(norad)
            continue
        sat = tles[norad]
        now = mean_elements_at(sat, epoch_jd)
        later = mean_elements_at(sat, epoch_jd + DECAY_FIT_DAYS)
        if now is None:
            report.sgp4_failed.append(norad)
            continue
        elements, n0 = now
        ndot = (later[1] - n0) / (DECAY_FIT_DAYS * 86400.0) if later is not None else 0.0
        rows.append(obj)
        els.append((*elements, ndot))

    table = table if table is not None else ObjectTable()
    if rows:
        cols = np.array(els, dtype=float).T
        table.add(
            ids=[o["norad_id"] for o in rows],
            names=[o["name"] for o in rows],
            category=[o.get("category", "active") for o in rows],
            mass_kg=[o.get("mass_kg", 100.0) for o in rows],
            radius_m=[o.get("radius_m", 1.0) for o in rows],
            elements=orbits.Elements(*cols),
            t_ref=0.0,
        )
    report.loaded = len(rows)
    return table, epoch_jd, report
