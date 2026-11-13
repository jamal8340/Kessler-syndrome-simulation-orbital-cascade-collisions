
### THIS IS KINDA TEMPORARY FILE. IT CAN BE USED TO CHECK THE FRAGMENTATION DIRECTLY FROM CLI. WILL PROBABLY DUMP IT WHEN THE LOGIC WILL BE INTEGRATED WITH UNITY

import json
import sys
from pathlib import Path
import os

CURRENT_FILE = Path(__file__).resolve()

PROJECT_ROOT = CURRENT_FILE.parents[2]

sys.path.insert(0, str(PROJECT_ROOT / "src"))

import fragmentation as fr  # noqa: E402

DEFAULT_INPUT = PROJECT_ROOT / "data" / "processed" / "close_approaches.json"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "processed" / "breakup_summary.json"

REQUIRED_FIELDS = (
    "mass_kg_1", "mass_kg_2", "radius_m_1", "radius_m_2",
    "velocity_kms_1", "velocity_kms_2",
)


def load_events(path):
    if not path.is_file():
        raise SystemExit(f"Brak pliku: {path}\nNajpierw uruchom close_approach.py.")
    events = json.loads(path.read_text(encoding="utf-8"))
    if events and any(f not in events[0] for f in REQUIRED_FIELDS):
        raise SystemExit(
            "close_approaches.json nie ma pól masy/promienia/prędkości.\n"
            "Uruchom ponownie close_approach.py po dodaniu pól z kroku 1."
        )
    return events


def select_events(events, pair, closest, event_type):
    if pair:
        wanted = {int(pair[0]), int(pair[1])}
        found = [e for e in events if {e["norad_id_1"], e["norad_id_2"]} == wanted]
        if not found:
            raise SystemExit(
                f"Pary {sorted(wanted)} nie ma na liście zbliżeń "
                "(nie przeszła progu odległości w analizowanym oknie)."
            )
        return [(found[0], found[0]["event_type"] != "PHYSICAL_COLLISION")]
    if closest:
        ranked = sorted(events, key=lambda e: e["miss_distance_km"])[:closest]
        return [(e, e["event_type"] != "PHYSICAL_COLLISION") for e in ranked]
    return [(e, False) for e in events if e["event_type"] == event_type]


def run_event(event, is_scenario, l_min_m, seed):
    a = fr.Body(event["mass_kg_1"], event["radius_m_1"], tuple(event["velocity_kms_1"]))
    b = fr.Body(event["mass_kg_2"], event["radius_m_2"], tuple(event["velocity_kms_2"]))
    r = fr.simulate_collision(a, b, l_min_m=l_min_m, seed=seed)

    largest = int(r.mass_kg.argmax()) if r.n_fragments else None
    return {
        "norad_id_1": event["norad_id_1"], "name_1": event["name_1"],
        "norad_id_2": event["norad_id_2"], "name_2": event["name_2"],
        "event_type": event["event_type"],
        "is_scenario": is_scenario,
        "tca_seconds_from_epoch": event["tca_seconds_from_epoch"],
        "miss_distance_km": event["miss_distance_km"],
        "relative_speed_kms": event["relative_speed_kms"],
        "encounter_altitude_km": event["encounter_altitude_km"],
        "mass_kg_1": a.mass_kg, "mass_kg_2": b.mass_kg,
        "catastrophic": bool(r.catastrophic),
        "emr_jg": round(r.emr_jg, 2),
        "n_expected": r.n_expected,
        "n_fragments": r.n_fragments,
        "tracked_mass_kg": round(float(r.mass_kg.sum()), 3),
        "untracked_mass_kg": round(r.remnant_mass_kg, 3),
        "largest_fragment_mass_kg": round(float(r.mass_kg[largest]), 3) if largest is not None else None,
        "largest_fragment_lc_m": round(float(r.lc_m[largest]), 3) if largest is not None else None,
        "dv_kinetic_fraction": float(f"{r.dv_kinetic_fraction:.3e}"),
        "l_min_m": l_min_m,
        "seed": seed,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--pair", nargs=2, metavar=("NORAD_1", "NORAD_2"))
    ap.add_argument("--closest", type=int, metavar="N")
    ap.add_argument("--event-type", default="PHYSICAL_COLLISION")
    ap.add_argument("--l-min", type=float, default=fr.L_MIN_DEFAULT_M)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    events = load_events(args.input)
    selected = select_events(events, args.pair, args.closest, args.event_type)

    if not selected:
        print(f"Brak zdarzeń typu {args.event_type} w {args.input.name} ({len(events)} zbliżeń).")
        print("Użyj --closest N albo --pair NORAD_1 NORAD_2, aby policzyć scenariusz.")
        return

    results = [
        run_event(e, scenario, args.l_min, args.seed + i)
        for i, (e, scenario) in enumerate(selected)
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=4, ensure_ascii=False), encoding="utf-8")

    print(f"=== FRAGMENTACJA: {len(results)} zdarzeń ===")
    for r in results:
        tag = "SCENARIUSZ" if r["is_scenario"] else r["event_type"]
        print(f"\n[{tag}] {r['name_1']} x {r['name_2']}  "
              f"(v_rel = {r['relative_speed_kms']} km/s, d_min = {r['miss_distance_km']} km)")
        print(f"  EMR = {r['emr_jg']} J/g -> {'katastroficzne' if r['catastrophic'] else 'niekatastroficzne'}")
        print(f"  Fragmenty >= {r['l_min_m']} m: {r['n_fragments']} (wzór: {r['n_expected']})")
        print(f"  Masa śledzona: {r['tracked_mass_kg']} kg, nieśledzona: {r['untracked_mass_kg']} kg")
    print(f"\nZapisano: {args.output}")


if __name__ == "__main__":
    main()