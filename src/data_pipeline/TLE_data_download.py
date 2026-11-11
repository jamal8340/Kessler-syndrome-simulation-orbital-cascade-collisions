import time
from pathlib import Path

import requests

CELESTRAK_GROUPS = {
    "active": "Aktywne satelity",
    "cosmos-2251-debris": "Odłamki Kosmos 2251 (kolizja z 2009 r.)",
    "iridium-33-debris": "Odłamki Iridium 33 (kolizja z 2009 r.)",
    "fengyun-1c-debris": "Odłamki Fengyun-1C (test ASAT 2007 r.)",
    "cosmos-1408-debris": "Odłamki Kosmos 1408 (test ASAT 2021 r.)",
}

BASE_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=tle"

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/122.0.0.0 Safari/537.36"
    )
}


def download_tle_groups():
    # Plik znajduje się w: <project_root>/src/data_pipeline/
    project_root = Path(__file__).resolve().parents[2]

    raw_dir = project_root / "data" / "raw" / "tle"
    raw_dir.mkdir(parents=True, exist_ok=True)

    legacy_file = project_root / "data" / "raw" / "active_satellites.txt"
    combined_entries = []

    for group_key, description in CELESTRAK_GROUPS.items():
        url = BASE_URL.format(group=group_key)
        group_file = raw_dir / f"{group_key}.tle"
        print(f"[INFO] Pobieranie grupy '{group_key}' ({description})...")

        try:
            response = requests.get(url, timeout=60, headers=HEADERS)
            response.raise_for_status()
            content = response.text.strip()

            if not content:
                print(f"[WARN] Grupa '{group_key}' zwróciła pustą odpowiedź.")
                continue

            with group_file.open("w", encoding="utf-8", newline="\n") as f:
                f.write(content + "\n")

            if group_key == "active":
                with legacy_file.open(
                    "w", encoding="utf-8", newline="\n"
                ) as f:
                    f.write(content + "\n")

            lines = [
                line.strip()
                for line in content.splitlines()
                if line.strip()
            ]
            num_objects = len(lines) // 3
            print(
                f"  -> Sukces! Pobrano ok. {num_objects} obiektów "
                f"(zapisano w {group_file})"
            )

            for i in range(0, len(lines) - 2, 3):
                if (
                    lines[i + 1].startswith("1 ")
                    and lines[i + 2].startswith("2 ")
                ):
                    combined_entries.append(
                        (group_key, lines[i], lines[i + 1], lines[i + 2])
                    )

            time.sleep(1.0)

        except requests.exceptions.RequestException as exc:
            print(
                f"[ERROR] Nie udało się pobrać grupy "
                f"'{group_key}': {exc}"
            )

    print(f"\nŁącznie pobrano {len(combined_entries)} wpisów TLE.")


if __name__ == "__main__":
    download_tle_groups()
