import json
import math
from datetime import datetime, timezone
from pathlib import Path

from sgp4.api import Satrec, jday

EARTH_RADIUS_KM = 6371.0
MIN_ALTITUDE_KM = 80.0
REFERENCE_FRAME = "TEME"

GROUP_FILES_ORDER = [
    ("cosmos-2251-debris", "cosmos-2251-debris.tle"),
    ("iridium-33-debris", "iridium-33-debris.tle"),
    ("fengyun-1c-debris", "fengyun-1c-debris.tle"),
    ("cosmos-1408-debris", "cosmos-1408-debris.tle"),
    ("active", "active.tle"),
]


def classify_object(name: str, source_group: str):
    upper_name = name.upper()
    if "R/B" in upper_name or "ROCKET" in upper_name:
        return "rocket_body", 2.5, 1500.0
    if "DEB" in upper_name or "debris" in source_group:
        return "debris", 0.15, 1.5
    return "active", 1.5, 500.0


def classify_orbit_regime(altitude_km: float, eccentricity: float) -> str:
    if eccentricity > 0.25:
        return "HEO"
    if altitude_km <= 2000.0:
        return "LEO"
    if 35000.0 <= altitude_km <= 36500.0:
        return "GEO"
    return "MEO"


def parse_tle_lines(lines, source_group, jd, fr, seen_norad_ids, stats):
    parsed_objects = []

    for i in range(len(lines) - 2):
        if lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
            stats["total_tle_records_read"] += 1
            name = lines[i]
            line1 = lines[i + 1]
            line2 = lines[i + 2]

            try:
                satellite = Satrec.twoline2rv(line1, line2)
                norad_id = int(satellite.satnum)

                if norad_id in seen_norad_ids:
                    stats["skipped_duplicate_norad"] += 1
                    continue

                error_code, r, v = satellite.sgp4(jd, fr)
                if error_code != 0:
                    stats["skipped_sgp4_error"] += 1
                    continue

                if any(not math.isfinite(val) for val in (*r, *v)):
                    stats["skipped_non_finite_state"] += 1
                    continue

                distance_from_center = math.sqrt(
                    r[0] ** 2 + r[1] ** 2 + r[2] ** 2
                )
                altitude_km = distance_from_center - EARTH_RADIUS_KM

                if altitude_km < MIN_ALTITUDE_KM:
                    stats["skipped_below_min_altitude"] += 1
                    continue

                speed_kms = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
                inclination_deg = math.degrees(satellite.inclo)
                eccentricity = float(satellite.ecco)
                bstar_drag = float(satellite.bstar)

                category, est_radius_m, est_mass_kg = classify_object(
                    name, source_group
                )
                orbit_regime = classify_orbit_regime(
                    altitude_km, eccentricity
                )

                seen_norad_ids.add(norad_id)

                parsed_objects.append({
                    "name": name,
                    "position": {
                        "x": round(r[0], 4),
                        "y": round(r[1], 4),
                        "z": round(r[2], 4),
                    },
                    "velocity": {
                        "vx": round(v[0], 5),
                        "vy": round(v[1], 5),
                        "vz": round(v[2], 5),
                    },
                    "norad_id": norad_id,
                    "category": category,
                    "source_group": source_group,
                    "orbit_regime": orbit_regime,
                    "altitude_km": round(altitude_km, 2),
                    "speed_kms": round(speed_kms, 4),
                    "inclination_deg": round(inclination_deg, 2),
                    "eccentricity": round(eccentricity, 6),
                    "bstar_drag": round(bstar_drag, 8),
                    "radius_m": est_radius_m,
                    "mass_kg": est_mass_kg,
                })
            except (ValueError, TypeError, OverflowError):
                stats["skipped_parse_exception"] += 1
                continue

    return parsed_objects


def process_tle_to_json():
    project_root = Path(__file__).resolve().parents[2]

    raw_dir = project_root / "data" / "raw" / "tle"
    legacy_input = project_root / "data" / "raw" / "active_satellites.txt"

    processed_dir = project_root / "data" / "processed"
    metadata_dir = project_root / "data" / "metadata"

    processed_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    output_file = processed_dir / "satellites_state.json"
    summary_file = metadata_dir / "simulation_summary.json"

    now = datetime.now(timezone.utc)
    jd, fr = jday(
        now.year, now.month, now.day,
        now.hour, now.minute, now.second
    )

    all_satellites = []
    seen_norad_ids = set()

    reconciliation_stats = {
        "total_tle_records_read": 0,
        "skipped_duplicate_norad": 0,
        "skipped_sgp4_error": 0,
        "skipped_non_finite_state": 0,
        "skipped_below_min_altitude": 0,
        "skipped_parse_exception": 0,
    }

    if raw_dir.is_dir():
        for group_name, filename in GROUP_FILES_ORDER:
            filepath = raw_dir / filename
            if filepath.is_file():
                with filepath.open("r", encoding="utf-8") as f:
                    lines = [
                        line.strip()
                        for line in f.read().splitlines()
                        if line.strip()
                    ]

                group_objs = parse_tle_lines(
                    lines, group_name, jd, fr, seen_norad_ids, reconciliation_stats
                )
                all_satellites.extend(group_objs)
                print(
                    f"Przetworzono grupę '{group_name}': "
                    f"dodano {len(group_objs)} unikalnych obiektów."
                )

    if not all_satellites and legacy_input.is_file():
        with legacy_input.open("r", encoding="utf-8") as f:
            lines = [
                line.strip()
                for line in f.read().splitlines()
                if line.strip()
            ]
        all_satellites = parse_tle_lines(
            lines, "active", jd, fr, seen_norad_ids, reconciliation_stats
        )

    if not all_satellites:
        print(
            "Nie przetworzono żadnych obiektów. "
            f"Sprawdź pliki TLE w: {raw_dir}"
        )
        return

    with output_file.open("w", encoding="utf-8") as f:
        json.dump(all_satellites, f, indent=4, ensure_ascii=False)

    total_skipped = (
        reconciliation_stats["skipped_duplicate_norad"]
        + reconciliation_stats["skipped_sgp4_error"]
        + reconciliation_stats["skipped_non_finite_state"]
        + reconciliation_stats["skipped_below_min_altitude"]
        + reconciliation_stats["skipped_parse_exception"]
    )
    reconciliation_stats["total_skipped"] = total_skipped
    reconciliation_stats["total_exported"] = len(all_satellites)

    summary = {
        "timestamp_utc": now.isoformat(),
        "reference_frame": REFERENCE_FRAME,
        "propagation_epoch_jd": round(jd + fr, 8),
        "minimum_altitude_threshold_km": MIN_ALTITUDE_KM,
        "total_objects": len(all_satellites),
        "record_reconciliation": reconciliation_stats,
        "by_category": {},
        "by_orbit_regime": {},
        "by_source_group": {},
        "leo_altitude_shells_km": {
            "200-400_km": 0,
            "400-600_km (Starlink_shell)": 0,
            "600-800_km": 0,
            "800-1000_km (Kessler_critical_zone)": 0,
            "1000-2000_km": 0,
        },
    }

    for obj in all_satellites:
        category = obj["category"]
        regime = obj["orbit_regime"]
        source_group = obj["source_group"]
        altitude = obj["altitude_km"]

        summary["by_category"][category] = (
            summary["by_category"].get(category, 0) + 1
        )
        summary["by_orbit_regime"][regime] = (
            summary["by_orbit_regime"].get(regime, 0) + 1
        )
        summary["by_source_group"][source_group] = (
            summary["by_source_group"].get(source_group, 0) + 1
        )

        if 200 <= altitude < 400:
            summary["leo_altitude_shells_km"]["200-400_km"] += 1
        elif 400 <= altitude < 600:
            summary["leo_altitude_shells_km"][
                "400-600_km (Starlink_shell)"
            ] += 1
        elif 600 <= altitude < 800:
            summary["leo_altitude_shells_km"]["600-800_km"] += 1
        elif 800 <= altitude < 1000:
            summary["leo_altitude_shells_km"][
                "800-1000_km (Kessler_critical_zone)"
            ] += 1
        elif 1000 <= altitude <= 2000:
            summary["leo_altitude_shells_km"]["1000-2000_km"] += 1

    with summary_file.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=4, ensure_ascii=False)

    print(f"\n=== PODSUMOWANIE PRZETWARZANIA SGP4 ({REFERENCE_FRAME}) ===")
    print(f"Wczytane rekordy TLE:     {reconciliation_stats['total_tle_records_read']}")
    print(f"Wyeksportowane obiekty:   {len(all_satellites)}")
    print(f"Odrzucone łącznie:        {total_skipped} (błąd SGP4/deorbitacja: {reconciliation_stats['skipped_sgp4_error']}, wys. < 80 km: {reconciliation_stats['skipped_below_min_altitude']})")
    print(f"Zapisano stany w:         {output_file}")
    print(f"Podsumowanie zapisano w:  {summary_file}")
    print("Rozkład kategorii:", summary["by_category"])
    print("Rozkład orbit:", summary["by_orbit_regime"])


if __name__ == "__main__":
    process_tle_to_json()
