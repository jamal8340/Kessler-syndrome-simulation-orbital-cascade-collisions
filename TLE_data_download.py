import requests
def download_active_satellites_tle():
    url = "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle"
    file_name = "active_satellites.txt"
    headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
        }
    try:
        response = requests.get(url, timeout=60, headers=headers)
        response.raise_for_status()

        with open(file_name, 'w', encoding="utf-8") as file:
            file.write(response.text)

        lines = response.text.strip().split("\n")
        if len(lines) >= 3:
            print("\ntest data: ")
            print('\n'.join(lines[:3]))

    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Error occured while downloading the satelite data.: {e}")

if __name__ == "__main__":
    download_active_satellites_tle()