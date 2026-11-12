
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt

# CONSTANTS


MU_EARTH = 398600.4418       # km^3/s^2
EARTH_RADIUS_KM = 6371.0

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT / "data" / "processed" / "satellites_state.json"
)

OUTPUT_FILE = (
    PROJECT_ROOT / "data" / "processed" / "orbital_elements.json"
)

FIGURES_DIR = PROJECT_ROOT / "results" / "figures"


# VECTOR OPERATIONS
def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    )


def magnitude(vector):
    return math.sqrt(dot(vector, vector))


def to_degrees(angle):
    return math.degrees(angle) % 360


# ORBITAL ELEMENTS

def calculate_orbital_elements(position, velocity):

    r = magnitude(position)
    v = magnitude(velocity)

    if r <= 0:
        raise ValueError("Invalid position vector")

    # Specific angular momentum
    h_vector = cross(position, velocity)
    h = magnitude(h_vector)

    if h <= 1e-12:
        raise ValueError("Angular momentum is zero")

    # Eccentricity vector
    vxh = cross(velocity, h_vector)

    eccentricity_vector = tuple(
        vxh[i] / MU_EARTH - position[i] / r
        for i in range(3)
    )

    eccentricity = magnitude(eccentricity_vector)

    # Specific orbital energy
    specific_energy = v**2 / 2 - MU_EARTH / r

    # Semi-major axis
    if specific_energy >= 0:
        raise ValueError("Object is not on a bound elliptical orbit")

    semi_major_axis = -MU_EARTH / (2 * specific_energy)

    # Inclination
    inclination = math.atan2(
        math.sqrt(h_vector[0]**2 + h_vector[1]**2),
        h_vector[2]
    )

    # Ascending node vector
    n_vector = (
        -h_vector[1],
        h_vector[0],
        0.0
    )

    n = magnitude(n_vector)

    # RAAN
    if n > 1e-10:
        raan = math.atan2(n_vector[1], n_vector[0]) % (2 * math.pi)
    else:
        raan = None

    # Argument of periapsis
    if n > 1e-10 and eccentricity > 1e-10:
        cos_argp = dot(n_vector, eccentricity_vector) / (
            n * eccentricity
        )

        sin_argp = dot(
            cross(n_vector, eccentricity_vector),
            h_vector
        ) / (n * eccentricity * h)

        arg_perigee = math.atan2(
            sin_argp, cos_argp
        ) % (2 * math.pi)
    else:
        arg_perigee = None

    # True anomaly
    if eccentricity > 1e-10:
        cos_nu = dot(eccentricity_vector, position) / (
            eccentricity * r
        )

        sin_nu = dot(
            cross(eccentricity_vector, position),
            h_vector
        ) / (eccentricity * r * h)

        true_anomaly = math.atan2(
            sin_nu, cos_nu
        ) % (2 * math.pi)
    else:
        true_anomaly = None

    # Perigee and apogee
    perigee_radius = semi_major_axis * (1 - eccentricity)
    apogee_radius = semi_major_axis * (1 + eccentricity)

    perigee_altitude = perigee_radius - EARTH_RADIUS_KM
    apogee_altitude = apogee_radius - EARTH_RADIUS_KM

    # Orbital period
    orbital_period = (
        2 * math.pi *
        math.sqrt(semi_major_axis**3 / MU_EARTH)
    ) / 60

    return {
        "semi_major_axis_km": semi_major_axis,
        "eccentricity": eccentricity,
        "inclination_deg": to_degrees(inclination),
        "raan_deg": to_degrees(raan) if raan is not None else None,
        "arg_perigee_deg": (
            to_degrees(arg_perigee)
            if arg_perigee is not None else None
        ),
        "true_anomaly_deg": (
            to_degrees(true_anomaly)
            if true_anomaly is not None else None
        ),
        "perigee_altitude_km": perigee_altitude,
        "apogee_altitude_km": apogee_altitude,
        "orbital_period_min": orbital_period,
        "specific_energy_km2_s2": specific_energy
    }



# PROCESS REAL POPULATION

def process_population():

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Nie znaleziono danych wejściowych: {INPUT_FILE}\n"
            "Najpierw uruchom process_TLE_to_json.py"
        )

    with INPUT_FILE.open("r", encoding="utf-8") as file:
        satellites = json.load(file)

    results = []

    skipped = 0

    print(f"Loaded objects: {len(satellites)}")

    for satellite in satellites:

        try:
            position = satellite["position"]
            velocity = satellite["velocity"]

            r = (
                position["x"],
                position["y"],
                position["z"]
            )

            v = (
                velocity["vx"],
                velocity["vy"],
                velocity["vz"]
            )

            elements = calculate_orbital_elements(r, v)

            result = {
                "norad_id": satellite["norad_id"],
                "name": satellite["name"],
                "category": satellite["category"],
                "source_group": satellite["source_group"],
                **elements
            }

            results.append(result)

        except (ValueError, KeyError, TypeError) as error:
            skipped += 1

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_FILE.open("w", encoding="utf-8") as file:
        json.dump(results, file, indent=4, ensure_ascii=False)

    print("\n=== ORBITAL MECHANICS ===")
    print(f"Input objects: {len(satellites)}")
    print(f"Successfully processed: {len(results)}")
    print(f"Skipped objects: {skipped}")
    print(f"Saved: {OUTPUT_FILE}")

    return results



# ORBITAL POPULATION PLOT


def generate_population_plot(results):

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Altitude of the semi-major axis
    altitudes = [
        obj["semi_major_axis_km"] - EARTH_RADIUS_KM
        for obj in results
    ]

    inclinations = [
        obj["inclination_deg"]
        for obj in results
    ]

    plt.figure(figsize=(13, 7))

    plt.scatter(
        altitudes,
        inclinations,
        s=3,
        alpha=0.35
    )

    plt.axvline(
        2000,
        linestyle="--",
        label="LEO upper boundary (2000 km)"
    )

    plt.axvline(
        35786,
        linestyle=":",
        label="Approximate GEO altitude"
    )

    plt.xlabel("Semi-major axis altitude [km]")
    plt.ylabel("Orbital inclination [degrees]")

    plt.title(
        "Orbital Population Distribution\n"
        f"{len(results)} catalogued objects"
    )

    plt.xlim(left=0)
    plt.ylim(0, 180)

    plt.grid(alpha=0.2)
    plt.legend()
    plt.tight_layout()

    output_plot = FIGURES_DIR / "orbital_population_map.png"

    plt.savefig(output_plot, dpi=200)
    plt.close()

    print(f"Plot saved: {output_plot}")



# MAIN


if __name__ == "__main__":

    orbital_results = process_population()

    generate_population_plot(orbital_results)

    print("\nOrbital analysis completed.")
