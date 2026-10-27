import json
import os
from datetime import datetime, timezone
from sgp4.api import Satrec, jday

def process_tle_to_json():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(script_dir, "active_satellites.txt")
    output_file = os.path.join(script_dir, "satellites_state.json")

    now = datetime.now(timezone.utc)
    jd, fr = jday(now.year, now.month, now.day, now.hour, now.minute, now.second)
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            # Używamy splitlines(), które radzi sobie ze znakami Windowsa, i usuwamy puste linijki
            lines = [line.strip() for line in f.read().splitlines() if line.strip()]
    except FileNotFoundError:
        print(f"[BŁĄD] Nie znaleziono pliku: {input_file}")
        return

    satellites = []
    
    # Przeszukujemy plik sprytniej - omijamy błędy przesunięcia
    for i in range(len(lines) - 2):
        # Sprawdzamy, czy dwie kolejne linijki to faktycznie poprawne dane orbity
        if lines[i+1].startswith('1 ') and lines[i+2].startswith('2 '):
            name = lines[i]
            line1 = lines[i+1]
            line2 = lines[i+2]
            
            try:
                satellite = Satrec.twoline2rv(line1, line2)
                error_code, r, v = satellite.sgp4(jd, fr) 
                
                # Tylko idealnie poprawne obliczenia trafiają do JSONa
                if error_code == 0:
                    satellites.append({
                        "name": name,
                        "position": {"x": r[0], "y": r[1], "z": r[2]},
                        "velocity": {"vx": v[0], "vy": v[1], "vz": v[2]}
                    })
            except Exception:
                continue

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(satellites, f, indent=4)
        
    print(f"[SUKCES] Przeliczono {len(satellites)} obiektów i zapisano w: {output_file}")

if __name__ == "__main__":
    process_tle_to_json()