import csv
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

from sgp4.api import Satrec

# Katalogi projektu:
# <project_root>/src/data_pipeline/validate_tle.py
PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = PROJECT_ROOT / "data" / "raw" / "tle"
OUTPUT_DIR = PROJECT_ROOT / "data" / "metadata"

STALE_AFTER_DAYS = 30


def tle_checksum_is_valid(line: str) -> bool:
    """Sprawdza standardową sumę kontrolną TLE (kolumny 1-69)."""
    if len(line) < 69 or not line[68].isdigit():
        return False

    checksum = 0
    for char in line[:68]:
        if char.isdigit():
            checksum += int(char)
        elif char == "-":
            checksum += 1

    return checksum % 10 == int(line[68])


def tle_epoch_datetime(satellite: Satrec) -> datetime:
    """Zwraca epokę TLE jako datę UTC."""
    year = int(satellite.epochyr)
    year += 2000 if year < 57 else 1900
    day_of_year = float(satellite.epochdays)

    epoch_start = datetime(year, 1, 1, tzinfo=timezone.utc)
    return epoch_start + timedelta(days=day_of_year - 1)


def extract_records(lines):
    """Odczytuje rekordy 3-liniowe (nazwa + linia 1 + linia 2)."""
    records = []
    malformed = 0
    i = 0

    while i < len(lines):
        current = lines[i].strip()

        if not current:
            i += 1
            continue

        if current.startswith("1 "):
            line1 = current
            name = "UNKNOWN"
            if i > 0 and not lines[i - 1].startswith(("1 ", "2 ")):
                name = lines[i - 1].strip() or "UNKNOWN"

            if i + 1 < len(lines) and lines[i + 1].strip().startswith("2 "):
                line2 = lines[i + 1].strip()
                records.append((name, line1, line2))
                i += 2
            else:
                malformed += 1
                i += 1
        else:
            i += 1

    return records, malformed


def validate_record(name, line1, line2, source_file, now):
    issues = []
    warnings = []
    norad_id = None
    epoch = None

    if len(line1) != 69:
        issues.append(f"linia 1 ma {len(line1)} znaków zamiast 69")
    if len(line2) != 69:
        issues.append(f"linia 2 ma {len(line2)} znaków zamiast 69")

    if not line1.startswith("1 "):
        issues.append("nieprawidłowy początek linii 1")
    if not line2.startswith("2 "):
        issues.append("nieprawidłowy początek linii 2")

    if len(line1) >= 7 and len(line2) >= 7:
        id1 = line1[2:7].strip()
        id2 = line2[2:7].strip()
        if id1.isdigit():
            norad_id = int(id1)
        if id1 != id2:
            issues.append(f"różne identyfikatory w liniach TLE: {id1} vs {id2}")

    if not tle_checksum_is_valid(line1):
        issues.append("niepoprawna suma kontrolna linii 1")
    if not tle_checksum_is_valid(line2):
        issues.append("niepoprawna suma kontrolna linii 2")

    satellite = None
    try:
        satellite = Satrec.twoline2rv(line1, line2)
        epoch = tle_epoch_datetime(satellite)

        age_days = (now - epoch).total_seconds() / 86400.0
        if age_days > STALE_AFTER_DAYS:
            warnings.append(f"TLE starsze niż {STALE_AFTER_DAYS} dni ({age_days:.1f} dni)")
        elif age_days < -1:
            warnings.append(f"epoka TLE jest w przyszłości ({-age_days:.1f} dni)")

        error_code, position, velocity = satellite.sgp4(
            satellite.jdsatepoch, satellite.jdsatepochF
        )
        if error_code != 0:
            issues.append(f"SGP4 zwrócił kod błędu {error_code}")
        elif not all(math.isfinite(value) for value in (*position, *velocity)):
            issues.append("SGP4 zwrócił niefinitywne wartości pozycji/prędkości")

    except Exception as exc:
        issues.append(f"nie można sparsować/propagować TLE: {type(exc).__name__}: {exc}")

    return {
        "source_file": source_file,
        "name": name,
        "norad_id": norad_id,
        "epoch_utc": epoch.isoformat() if epoch else None,
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
    }


def main():
    now = datetime.now(timezone.utc)

    if not INPUT_DIR.is_dir():
        raise SystemExit(f"Nie znaleziono katalogu z TLE: {INPUT_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []
    duplicate_ids = {}
    total_records = 0
    malformed_records = 0

    tle_files = sorted(INPUT_DIR.glob("*.tle"))
    if not tle_files:
        raise SystemExit(f"Nie znaleziono plików .tle w: {INPUT_DIR}")

    for tle_file in tle_files:
        lines = [
            line.strip()
            for line in tle_file.read_text(
                encoding="utf-8", errors="replace"
            ).splitlines()
            if line.strip()
        ]

        records, malformed = extract_records(lines)
        malformed_records += malformed
        print(f"{tle_file.name}: znaleziono {len(records)} rekordów; niekompletnych: {malformed}")

        for name, line1, line2 in records:
            result = validate_record(name, line1, line2, tle_file.name, now)
            all_results.append(result)
            total_records += 1

            if result["norad_id"] is not None:
                duplicate_ids.setdefault(result["norad_id"], []).append(tle_file.name)

    valid_count = sum(1 for item in all_results if item["valid"])
    invalid_count = total_records - valid_count
    warning_count = sum(1 for item in all_results if item["warnings"])

    duplicate_list = [
        {"norad_id": norad_id, "files": files}
        for norad_id, files in duplicate_ids.items()
        if len(files) > 1
    ]

    report = {
        "generated_utc": now.isoformat(),
        "input_directory": str(INPUT_DIR),
        "files_checked": len(tle_files),
        "records_found": total_records,
        "valid_records": valid_count,
        "invalid_records": invalid_count,
        "records_with_warnings": warning_count,
        "incomplete_records": malformed_records,
        "duplicate_norad_ids_across_files": duplicate_list,
        "results": all_results,
    }

    json_path = OUTPUT_DIR / "tle_validation_report.json"
    csv_path = OUTPUT_DIR / "tle_validation_details.csv"

    json_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    with csv_path.open("w", newline="", encoding="utf-8-sig") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=[
                "source_file", "name", "norad_id", "epoch_utc",
                "valid", "issues", "warnings",
            ],
        )
        writer.writeheader()
        for item in all_results:
            row = item.copy()
            row["issues"] = "; ".join(row["issues"])
            row["warnings"] = "; ".join(row["warnings"])
            writer.writerow(row)

    print("\n=== PODSUMOWANIE WALIDACJI TLE ===")
    print(f"Pliki sprawdzone:       {len(tle_files)}")
    print(f"Rekordy znalezione:     {total_records}")
    print(f"Poprawne rekordy:       {valid_count}")
    print(f"Niepoprawne rekordy:    {invalid_count}")
    print(f"Rekordy z ostrzeżeniem: {warning_count}")
    print(f"Niekompletne rekordy:   {malformed_records}")
    print(f"Powtórzone NORAD ID:    {len(duplicate_list)}")
    print(f"\nRaport JSON: {json_path}")
    print(f"Raport CSV:  {csv_path}")


if __name__ == "__main__":
    main()
