import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt

# CONSTANTS

MU_EARTH = 398600.4418       # km^3/s^2
EARTH_RADIUS_KM = 6371.0

DEFAULT_SCREENING_THRESHOLD_KM = 10.0
CRITICAL_CONJUNCTION_KM = 2.0
CELL_SIZE_KM = 80.0

PROPAGATION_STEP_SECONDS = 30.0
PROPAGATION_STEPS_COUNT = 20

MIN_RELATIVE_SPEED_KMS = 0.5
CATASTROPHIC_ENERGY_RATIO_JG = 40.0  

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "data" / "processed" / "satellites_state.json"
OUTPUT_JSON = PROJECT_ROOT / "data" / "processed" / "close_approaches.json"
REPORT_JSON = PROJECT_ROOT / "data" / "metadata" / "collision_detection_report.json"
REPORT_CSV = PROJECT_ROOT / "data" / "metadata" / "close_approaches_table.csv"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"

FORWARD_NEIGHBOUR_OFFSETS = [
    (dx, dy, dz)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    for dz in (-1, 0, 1)
    if (dx > 0) or (dx == 0 and dy > 0) or (dx == 0 and dy == 0 and dz > 0)
]



def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add_scaled(a, b, scale):
    return (
        a[0] + b[0] * scale,
        a[1] + b[1] * scale,
        a[2] + b[2] * scale,
    )


def magnitude(v):
    return math.sqrt(dot(v, v))


def compute_tca_and_miss_distance(r_a, v_a, r_b, v_b, dt_window):
    dr0 = sub(r_a, r_b)
    dv = sub(v_a, v_b)

    rel_speed_sq = dot(dv, dv)
    instantaneous_dist = magnitude(dr0)

    if rel_speed_sq <= 1e-12:
        return 0.0, instantaneous_dist, 0.0, instantaneous_dist

    rel_speed = math.sqrt(rel_speed_sq)

    t_unconstrained = -dot(dr0, dv) / rel_speed_sq
    t_tca = max(0.0, min(dt_window, t_unconstrained))

    dr_at_tca = add_scaled(dr0, dv, t_tca)
    miss_distance_km = magnitude(dr_at_tca)

    return t_tca, miss_distance_km, rel_speed, instantaneous_dist


def compute_collision_energetics(mass_a_kg, mass_b_kg, rel_speed_kms):
    if mass_a_kg <= 0 or mass_b_kg <= 0:
        return 0.0, 0.0, False

    rel_speed_ms = rel_speed_kms * 1000.0
    reduced_mass = (mass_a_kg * mass_b_kg) / (mass_a_kg + mass_b_kg)

    kinetic_energy_joules = 0.5 * reduced_mass * (rel_speed_ms ** 2)
    kinetic_energy_mj = kinetic_energy_joules / 1e6

    target_mass_kg = max(mass_a_kg, mass_b_kg)
    emr_joules_per_gram = kinetic_energy_joules / (target_mass_kg * 1000.0)
    is_catastrophic = emr_joules_per_gram >= CATASTROPHIC_ENERGY_RATIO_JG

    return kinetic_energy_mj, emr_joules_per_gram, is_catastrophic


def is_same_constellation_formation(obj_a, obj_b):
    if obj_a["category"] != "active" or obj_b["category"] != "active":
        return False

    name_a = obj_a["name"].upper()
    name_b = obj_b["name"].upper()

    for prefix in ("STARLINK", "ONEWEB", "IRIDIUM", "FLOCK"):
        if name_a.startswith(prefix) and name_b.startswith(prefix):
            return True

    return False



def build_spatial_grid(objects, cell_size_km):
    grid = {}
    for idx, obj in enumerate(objects):
        r = obj["r"]
        cell = (
            math.floor(r[0] / cell_size_km),
            math.floor(r[1] / cell_size_km),
            math.floor(r[2] / cell_size_km),
        )
        if cell not in grid:
            grid[cell] = [idx]
        else:
            grid[cell].append(idx)
    return grid


def screen_candidate_pairs(objects, cell_size_km):
    grid = build_spatial_grid(objects, cell_size_km)
    candidate_pairs = []

    for cell, indices in grid.items():
        n = len(indices)
        for a_idx in range(n):
            i = indices[a_idx]
            for b_idx in range(a_idx + 1, n):
                j = indices[b_idx]
                candidate_pairs.append((i, j))

        cx, cy, cz = cell
        for dx, dy, dz in FORWARD_NEIGHBOUR_OFFSETS:
            neighbour_cell = (cx + dx, cy + dy, cz + dz)
            neighbour_indices = grid.get(neighbour_cell)
            if neighbour_indices:
                for i in indices:
                    for j in neighbour_indices:
                        candidate_pairs.append((i, j))

    return candidate_pairs



def propagate_step_rk2(objects, dt):
    min_radius_sq = (EARTH_RADIUS_KM + 80.0) ** 2

    for obj in objects:
        r0 = obj["r"]
        v0 = obj["v"]

        r0_sq = dot(r0, r0)
        if r0_sq < min_radius_sq:
            continue

        r0_mag = math.sqrt(r0_sq)
        factor0 = -MU_EARTH / (r0_sq * r0_mag)
        a0 = (r0[0] * factor0, r0[1] * factor0, r0[2] * factor0)

        r_mid = add_scaled(r0, v0, 0.5 * dt)
        v_mid = add_scaled(v0, a0, 0.5 * dt)

        r_mid_sq = dot(r_mid, r_mid)
        r_mid_mag = math.sqrt(r_mid_sq)
        factor_mid = -MU_EARTH / (r_mid_sq * r_mid_mag)
        a_mid = (
            r_mid[0] * factor_mid,
            r_mid[1] * factor_mid,
            r_mid[2] * factor_mid,
        )

        obj["r"] = add_scaled(r0, v_mid, dt)
        obj["v"] = add_scaled(v0, a_mid, dt)



def run_validation_tests():
    print("=== URUCHAMIANIE TESTÓW WALIDACYJNYCH DETEKCJI KOLIZJI ===")

    r_a = (-10.0, 0.0, 7000.0)
    v_a = (1.0, 0.0, 0.0)
    r_b = (10.0, 0.0, 7000.0)
    v_b = (-1.0, 0.0, 0.0)
    t_tca, d_min, v_rel, d_inst = compute_tca_and_miss_distance(
        r_a, v_a, r_b, v_b, dt_window=20.0
    )
    assert abs(t_tca - 10.0) < 1e-6, f"Test 1 TCA error: {t_tca}"
    assert abs(d_min - 0.0) < 1e-6, f"Test 1 d_min error: {d_min}"
    assert abs(v_rel - 2.0) < 1e-6, f"Test 1 v_rel error: {v_rel}"
    print(
        f"[PASS] Test 1 (Zderzenie czołowe): d_start={d_inst:.1f} km -> "
        f"TCA={t_tca:.2f} s, d_min={d_min:.4f} km, v_rel={v_rel:.2f} km/s"
    )

    r_a2 = (-60.0, 1.2, 7100.0)
    v_a2 = (6.0, 0.0, 0.0)
    r_b2 = (60.0, 0.0, 7100.0)
    v_b2 = (-6.0, 0.0, 0.0)
    t_tca2, d_min2, v_rel2, d_inst2 = compute_tca_and_miss_distance(
        r_a2, v_a2, r_b2, v_b2, dt_window=20.0
    )
    assert abs(t_tca2 - 10.0) < 1e-6, f"Test 2 TCA error: {t_tca2}"
    assert abs(d_min2 - 1.2) < 1e-6, f"Test 2 d_min error: {d_min2}"
    print(
        f"[PASS] Test 2 (Koniunkcja wewnątrz kroku): d_start={d_inst2:.2f} km -> "
        f"TCA={t_tca2:.2f} s, d_min={d_min2:.4f} km (wykryto mimo dużego kroku 20 s)"
    )

    mock_objects = [
        {"r": (CELL_SIZE_KM - 0.1, 0.0, 7000.0)},
        {"r": (CELL_SIZE_KM + 0.1, 0.0, 7000.0)},
    ]
    pairs = screen_candidate_pairs(mock_objects, CELL_SIZE_KM)
    assert len(pairs) == 1, "Test 3 failed: boundary pair missed by spatial grid!"
    print(
        f"[PASS] Test 3 (Granica komórek 3D): obiekty na {CELL_SIZE_KM - 0.1} km "
        f"oraz {CELL_SIZE_KM + 0.1} km poprawnie sparowane przez sąsiedztwo 27 komórek.\n"
    )



def run_catalogue_screening():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku wejściowego: {INPUT_FILE}\n"
            "Najpierw uruchom src/data_pipeline/process_TLE_to_json.py"
        )

    with INPUT_FILE.open("r", encoding="utf-8") as f:
        raw_satellites = json.load(f)

    objects = []
    for sat in raw_satellites:
        pos = sat["position"]
        vel = sat["velocity"]
        objects.append({
            "norad_id": sat["norad_id"],
            "name": sat["name"],
            "category": sat["category"],
            "orbit_regime": sat.get("orbit_regime", "LEO"),
            "radius_m": sat.get("radius_m", 1.0),
            "mass_kg": sat.get("mass_kg", 100.0),
            "r": (float(pos["x"]), float(pos["y"]), float(pos["z"])),
            "v": (float(vel["vx"]), float(vel["vy"]), float(vel["vz"])),
        })

    print(f"Wczytano obiektów do analizy koniunkcji: {len(objects)}")
    print(
        f"Horyzont analizy: {PROPAGATION_STEPS_COUNT * PROPAGATION_STEP_SECONDS:.0f} s "
        f"({PROPAGATION_STEPS_COUNT} kroków po {PROPAGATION_STEP_SECONDS:.0f} s, "
        f"siatka 3D: {CELL_SIZE_KM:.0f} km, próg: {DEFAULT_SCREENING_THRESHOLD_KM:.1f} km)"
    )

    best_encounters_by_pair = {}
    total_candidate_pairs_checked = 0
    ignored_constellation_pairs = 0

    for step_idx in range(PROPAGATION_STEPS_COUNT):
        step_start_time_s = step_idx * PROPAGATION_STEP_SECONDS
        candidate_pairs = screen_candidate_pairs(objects, CELL_SIZE_KM)
        total_candidate_pairs_checked += len(candidate_pairs)

        for i, j in candidate_pairs:
            obj_a = objects[i]
            obj_b = objects[j]

            if is_same_constellation_formation(obj_a, obj_b):
                ignored_constellation_pairs += 1
                continue

            t_local, d_min, rel_speed, d_inst = compute_tca_and_miss_distance(
                obj_a["r"], obj_a["v"],
                obj_b["r"], obj_b["v"],
                PROPAGATION_STEP_SECONDS
            )

            if d_min <= DEFAULT_SCREENING_THRESHOLD_KM and rel_speed >= MIN_RELATIVE_SPEED_KMS:
                pair_key = (
                    min(obj_a["norad_id"], obj_b["norad_id"]),
                    max(obj_a["norad_id"], obj_b["norad_id"]),
                )

                existing = best_encounters_by_pair.get(pair_key)
                if existing is None or d_min < existing["miss_distance_km"]:
                    r_a_tca = add_scaled(obj_a["r"], obj_a["v"], t_local)
                    r_b_tca = add_scaled(obj_b["r"], obj_b["v"], t_local)
                    midpoint = (
                        0.5 * (r_a_tca[0] + r_b_tca[0]),
                        0.5 * (r_a_tca[1] + r_b_tca[1]),
                        0.5 * (r_a_tca[2] + r_b_tca[2]),
                    )
                    altitude_km = magnitude(midpoint) - EARTH_RADIUS_KM

                    hard_body_radius_km = (
                        obj_a["radius_m"] + obj_b["radius_m"]
                    ) / 1000.0

                    if d_min <= hard_body_radius_km:
                        event_type = "PHYSICAL_COLLISION"
                    elif d_min <= CRITICAL_CONJUNCTION_KM:
                        event_type = "CRITICAL_CONJUNCTION"
                    else:
                        event_type = "CLOSE_APPROACH"

                    ke_mj, emr_jg, catastrophic = compute_collision_energetics(
                        obj_a["mass_kg"], obj_b["mass_kg"], rel_speed
                    )

                    best_encounters_by_pair[pair_key] = {
                        "norad_id_1": obj_a["norad_id"],
                        "name_1": obj_a["name"],
                        "category_1": obj_a["category"],
                        "norad_id_2": obj_b["norad_id"],
                        "name_2": obj_b["name"],
                        "category_2": obj_b["category"],
                        "event_type": event_type,
                        "tca_seconds_from_epoch": round(step_start_time_s + t_local, 3),
                        "miss_distance_km": round(d_min, 4),
                        "initial_step_distance_km": round(d_inst, 4),
                        "hard_body_radius_km": round(hard_body_radius_km, 6),
                        "relative_speed_kms": round(rel_speed, 4),
                        "encounter_altitude_km": round(altitude_km, 2),
                        "encounter_position_km": {
                            "x": round(midpoint[0], 3),
                            "y": round(midpoint[1], 3),
                            "z": round(midpoint[2], 3),
                        },
                        "potential_kinetic_energy_mj": round(ke_mj, 2),
                        "energy_to_mass_ratio_jg": round(emr_jg, 2),
                        "would_be_catastrophic": catastrophic,
                        "debris_involved": (
                            obj_a["category"] != "active"
                            or obj_b["category"] != "active"
                        ),
                    }

        if step_idx < PROPAGATION_STEPS_COUNT - 1:
            propagate_step_rk2(objects, PROPAGATION_STEP_SECONDS)

    encounters = sorted(
        best_encounters_by_pair.values(),
        key=lambda x: x["miss_distance_km"]
    )

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(encounters, f, indent=4, ensure_ascii=False)

    if encounters:
        csv_fields = [
            "norad_id_1", "name_1", "category_1",
            "norad_id_2", "name_2", "category_2",
            "event_type", "tca_seconds_from_epoch",
            "miss_distance_km", "relative_speed_kms",
            "encounter_altitude_km", "potential_kinetic_energy_mj",
            "energy_to_mass_ratio_jg", "would_be_catastrophic", "debris_involved",
        ]
        with REPORT_CSV.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(encounters)

    physical_collisions = sum(
        1 for e in encounters if e["event_type"] == "PHYSICAL_COLLISION"
    )
    critical_conjunctions = sum(
        1 for e in encounters if e["event_type"] == "CRITICAL_CONJUNCTION"
    )
    debris_encounters = sum(1 for e in encounters if e["debris_involved"])
    catastrophic_potential = sum(
        1 for e in encounters if e["would_be_catastrophic"]
    )

    report = {
        "input_objects_count": len(objects),
        "simulation_horizon_seconds": PROPAGATION_STEPS_COUNT * PROPAGATION_STEP_SECONDS,
        "time_step_seconds": PROPAGATION_STEP_SECONDS,
        "spatial_cell_size_km": CELL_SIZE_KM,
        "screening_threshold_km": DEFAULT_SCREENING_THRESHOLD_KM,
        "critical_conjunction_threshold_km": CRITICAL_CONJUNCTION_KM,
        "total_candidate_pairs_evaluated": total_candidate_pairs_checked,
        "ignored_same_constellation_pairs": ignored_constellation_pairs,
        "unique_close_approaches_found": len(encounters),
        "physical_collisions_count": physical_collisions,
        "critical_conjunctions_under_2km": critical_conjunctions,
        "encounters_involving_debris": debris_encounters,
        "catastrophic_potential_count": catastrophic_potential,
        "closest_5_encounters": encounters[:5],
    }

    with REPORT_JSON.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=4, ensure_ascii=False)

    print("=== WYNIKI DETEKCJI BLISKICH PODEJŚĆ (STAGE 5) ===")
    print(f"Przebadane pary w siatce 27-komórkowej: {total_candidate_pairs_checked}")
    print(f"Wykryte unikalne zbliżenia (<= {DEFAULT_SCREENING_THRESHOLD_KM} km): {len(encounters)}")
    print(f"Krytyczne koniunkcje (<= {CRITICAL_CONJUNCTION_KM} km):        {critical_conjunctions}")
    print(f"Bezpośrednie kolizje fizyczne:          {physical_collisions}")
    print(f"Zbliżenia z udziałem śmieci (debris):   {debris_encounters}")
    print(f"Potencjalnie katastroficzne (>40 J/g):  {catastrophic_potential}")
    if encounters:
        top = encounters[0]
        print(
            f"\nNajbliższe przelot: {top['name_1']} x {top['name_2']} -> "
            f"d_min = {top['miss_distance_km']} km (TCA = {top['tca_seconds_from_epoch']} s, "
            f"v_rel = {top['relative_speed_kms']} km/s, wys. = {top['encounter_altitude_km']} km)"
        )
    print(f"\nZapisano pełną listę w: {OUTPUT_JSON}")
    print(f"Zapisano raport JSON w: {REPORT_JSON}")
    print(f"Zapisano tabelę CSV w:  {REPORT_CSV}")

    return encounters



def generate_conjunctions_plot(encounters):
    if not encounters:
        return

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    debris_events = [e for e in encounters if e["debris_involved"]]
    active_only_events = [e for e in encounters if not e["debris_involved"]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    if active_only_events:
        ax1.scatter(
            [e["encounter_altitude_km"] for e in active_only_events],
            [e["miss_distance_km"] for e in active_only_events],
            s=18, alpha=0.6, label="Active vs Active", color="#1f77b4"
        )
    if debris_events:
        ax1.scatter(
            [e["encounter_altitude_km"] for e in debris_events],
            [e["miss_distance_km"] for e in debris_events],
            s=22, alpha=0.75, label="Debris / Rocket Body involved", color="#d62728"
        )

    ax1.axhline(
        CRITICAL_CONJUNCTION_KM,
        color="red", linestyle="--",
        label=f"Critical threshold ({CRITICAL_CONJUNCTION_KM} km)"
    )
    ax1.set_xlabel("Encounter Altitude [km]")
    ax1.set_ylabel("Minimum Miss Distance at TCA [km]")
    ax1.set_title("Close Approaches: Altitude vs. Miss Distance")
    ax1.grid(alpha=0.25)
    ax1.legend()

    ax2.scatter(
        [e["encounter_altitude_km"] for e in encounters],
        [e["relative_speed_kms"] for e in encounters],
        c=[e["miss_distance_km"] for e in encounters],
        cmap="viridis_r", s=20, alpha=0.75
    )
    ax2.set_xlabel("Encounter Altitude [km]")
    ax2.set_ylabel("Relative Impact Speed [km/s]")
    ax2.set_title("Relative Encounter Velocity across Orbital Altitudes")
    ax2.grid(alpha=0.25)

    plt.tight_layout()
    plot_path = FIGURES_DIR / "close_approaches_analysis.png"
    plt.savefig(plot_path, dpi=200)
    plt.close()

    print(f"Wykres analizy koniunkcji zapisano w: {plot_path}")


if __name__ == "__main__":
    run_validation_tests()
    detected_encounters = run_catalogue_screening()
    generate_conjunctions_plot(detected_encounters)