"""Kopiuje wspólny zbiór danych (data/processed/satellites_state.json) do projektu Unity.

Obok pliku stanów zapisuje dataset_info.json (epoka, układ odniesienia, liczba obiektów,
SHA-256), żeby obie strony mogły sprawdzić, że używają tego samego zbioru.
Istniejący plik w Unity, jeśli się różni, jest kopiowany do <unity>/DatasetBackups/ (poza Assets).
"""

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
STATE_FILE = PROJECT_ROOT / "data" / "processed" / "satellites_state.json"
SUMMARY_FILE = PROJECT_ROOT / "data" / "metadata" / "simulation_summary.json"
DEFAULT_UNITY_PROJECT = PROJECT_ROOT.parent / "KesslerSimulation"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def export(unity_project):
    target_dir = unity_project / "Assets" / "StreamingAssets"
    if not (unity_project / "Assets").is_dir():
        raise SystemExit(f"To nie jest projekt Unity (brak katalogu Assets): {unity_project}")
    target_dir.mkdir(parents=True, exist_ok=True)

    summary = json.loads(SUMMARY_FILE.read_text(encoding="utf-8"))
    count = len(json.loads(STATE_FILE.read_text(encoding="utf-8")))
    if count != summary["total_objects"]:
        raise SystemExit(
            f"Niespójne dane: {STATE_FILE.name} ma {count} obiektów, "
            f"a {SUMMARY_FILE.name} podaje {summary['total_objects']}. Uruchom ponownie process_TLE_to_json.py."
        )

    digest = sha256(STATE_FILE)
    target = target_dir / "satellites_state.json"

    if target.is_file() and sha256(target) != digest:
        backup_dir = unity_project / "DatasetBackups"
        backup_dir.mkdir(exist_ok=True)
        stamp = datetime.fromtimestamp(target.stat().st_mtime, timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        backup = backup_dir / f"satellites_state_{stamp}.json"
        shutil.copy2(target, backup)
        print(f"Kopia zapasowa poprzedniego pliku Unity: {backup}")

    shutil.copyfile(STATE_FILE, target)

    info = {
        "epoch_utc": summary["timestamp_utc"],
        "epoch_jd": summary["propagation_epoch_jd"],
        "reference_frame": summary["reference_frame"],
        "units": {"position": "km", "velocity": "km/s"},
        "object_count": count,
        "sha256": digest,
        "exported_utc": datetime.now(timezone.utc).isoformat(),
    }
    (target_dir / "dataset_info.json").write_text(json.dumps(info, indent=4), encoding="utf-8")

    print(f"Wyeksportowano {count} obiektów (epoka {info['epoch_utc']}, {info['reference_frame']}) do: {target}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--unity-project", type=Path, default=DEFAULT_UNITY_PROJECT)
    export(parser.parse_args().unity_project.resolve())
