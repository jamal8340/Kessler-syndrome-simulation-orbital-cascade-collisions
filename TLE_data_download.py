import os
import time
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
    script_dir = os.path.dirname(os.path.abspath(__file__))
    raw_dir = os.path.join(script_dir, "tle_raw")
    os.makedirs(raw_dir, exist_ok=True)

    combined_entries = []

    for group_key, description in CELESTRAK_GROUPS.items():
        url = BASE_URL.format(group=group_key)
        group_file = os.path.join(raw_dir, f"{group_key}.tle")
        print(f"[INFO] Pobieranie grupy '{group_key}' ({description})...")

        try:
            response = requests.get(url, timeout=60, headers=HEADERS)
            response.raise_for_status()
            content = response.text.strip()

            with open(group_file, "w", encoding="utf-8") as f:
                f.write(content + "\n")

            if group_key == "active":
                legacy_file = os.path.join(script_dir, "active_satellites.txt")
                with open(legacy_file, "w", encoding="utf-8") as f:
                    f.write(content + "\n")

            lines = [line.strip() for line in content.splitlines() if line.strip()]
            num_objects = len(lines) // 3
            print(f"  -> Sukces! Pobrano ok. {num_objects} obiektów (zapisano w {group_file})")

            for i in range(0, len(lines) - 2, 3):
                if lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
                    combined_entries.append((group_key, lines[i], lines[i + 1], lines[i + 2]))

            time.sleep(1.0)

        except requests.exceptions.RequestException as e:
            print(f"Nie udało się pobrać grupy '{group_key}': {e}")

    print(f"\nŁącznie pobrano {len(combined_entries)} wpisów TLE.")


if __name__ == "__main__":
    download_tle_groups()