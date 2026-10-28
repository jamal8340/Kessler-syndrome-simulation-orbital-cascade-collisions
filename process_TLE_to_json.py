import json
import math
import os
from datetime import datetime, timezone
from sgp4.api import Satrec, jday

EARTH_RADIUS_KM = 6371.0

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


def parse_tle_lines(lines, source_group, jd, fr, seen_norad_ids):
    parsed_objects = []

    for i in range(len(lines) - 2):
        if lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
            name = lines[i]
            line1 = lines[i + 1]
            line2 = lines[i + 2]

            try:
                satellite = Satrec.twoline2rv(line1, line2)
                norad_id = int(satellite.satnum)

                if norad_id in seen_norad_ids:
                    continue

                error_code, r, v = satellite.sgp4(jd, fr)
                if error_code != 0:
                    continue

                if any(math.isnan(val) for val in r) or any(math.isnan(val) for val in v):
                    continue

                distance_from_center = math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
                altitude_km = distance_from_center - EARTH_RADIUS_KM

                if altitude_km < 80.0:
                    continue

                speed_kms = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
                inclination_deg = math.degrees(satellite.inclo)
                eccentricity = float(satellite.ecco)
                bstar_drag = float(satellite.bstar)

                category, est_radius_m, est_mass_kg = classify_object(name, source_group)
                orbit_regime = classify_orbit_regime(altitude_km, eccentricity)

                seen_norad_ids.add(norad_id)

                parsed_objects.append({
                    "name": name,
                    "position": {"x": round(r[0], 4), "y": round(r[1], 4), "z": round(r[2], 4)},
                    "velocity": {"vx": round(v[0], 5), "vy": round(v[1], 5), "vz": round(v[2], 5)},
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
            except Exception:
                continue

    return parsed_objects


def process_tle_to_json():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    raw_dir = os.path.join(script_dir, "tle_raw")
    output_file = os.path.join(script_dir, "satellites_state.json")
    summary_file = os.path.join(script_dir, "simulation_summary.json")

    now = datetime.now(timezone.utc)
    jd, fr = jday(now.year, now.month, now.day, now.hour, now.minute, now.second)

    all_satellites = []
    seen_norad_ids = set()

    if os.path.isdir(raw_dir):
        for group_name, filename in GROUP_FILES_ORDER:
            filepath = os.path.join(raw_dir, filename)
            if os.path.exists(filepath):
                with open(filepath, "r", encoding="utf-8") as f:
                    lines = [line.strip() for line in f.read().splitlines() if line.strip()]
                group_objs = parse_tle_lines(lines, group_name, jd, fr, seen_norad_ids)
                all_satellites.extend(group_objs)
                print(f"Przetworzono grupę '{group_name}': dodano {len(group_objs)} unikalnych obiektów.")

    if not all_satellites:
        legacy_input = os.path.join(script_dir, "active_satellites.txt")
        if os.path.exists(legacy_input):
            with open(legacy_input, "r", encoding="utf-8") as f:
                lines = [line.strip() for line in f.read().splitlines() if line.strip()]
            all_satellites = parse_tle_lines(lines, "active", jd, fr, seen_norad_ids)
        else:
            print("Brak plików TLE! Uruchom najpierw TLE_data_download.py.")
            return

    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(all_satellites, f, indent=4)

    summary = {
        "timestamp_utc": now.isoformat(),
        "total_objects": len(all_satellites),
        "by_category": {},
        "by_orbit_regime": {},
        "by_source_group": {},
        "leo_altitude_ shells_km": {
            "200-400_km": 0,
            "400-600_km (Starlink_shell)": 0,
            "600-800_km": 0,
            "800-1000_km (Kessler_critical_zone)": 0,
            "1000-2000_km": 0,
        },
    }

    for obj in all_satellites:
        cat = obj["category"]
        reg = obj["orbit_regime"]
        grp = obj["source_group"]
        alt = obj["altitude_km"]

        summary["by_category"][cat] = summary["by_category"].get(cat, 0) + 1
        summary["by_orbit_regime"][reg] = summary["by_orbit_regime"].get(reg, 0) + 1
        summary["by_source_group"][grp] = summary["by_source_group"].get(grp, 0) + 1

        if 200 <= alt < 400:
            summary["leo_altitude_ shells_km"]["200-400_km"] += 1
        elif 400 <= alt < 600:
            summary["leo_altitude_ shells_km"]["400-600_km (Starlink_shell)"] += 1
        elif 600 <= alt < 800:
            summary["leo_altitude_ shells_km"]["600-800_km"] += 1
        elif 800 <= alt < 1000:
            summary["leo_altitude_ shells_km"]["800-1000_km (Kessler_critical_zone)"] += 1
        elif 1000 <= alt <= 2000:
            summary["leo_altitude_ shells_km"]["1000-2000_km"] += 1

    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=4)

    print(f"\nZapisano {len(all_satellites)} obiektów w: {output_file}")
    print(f"Podsumowanie statystyczne zapisano w: {summary_file}")
    print("Rozkład kategorii:", summary["by_category"])
    print("Rozkład orbit:", summary["by_orbit_regime"])


if __name__ == "__main__":
    process_tle_to_json()