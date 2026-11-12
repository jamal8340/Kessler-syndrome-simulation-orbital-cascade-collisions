# Kessler Syndrome Simulation

> **Status:** Research prototype under development\
> **Main technologies:** Python, SGP4, Matplotlib, Unity 6\
> **Goal:** orbital population analysis, close-approach detection,
> collision modelling and cascading debris evolution.

## 0. What's New in the Current Development Stage

Recent work: - reorganized the Python project directories; - corrected
paths in the TLE download and processing pipeline; - added TLE
validation with JSON/CSV reports; - fixed parsing of CelesTrak files
containing blank lines; - implemented osculating orbital element
calculations; - processed 18,479 objects with zero skipped in the latest
reported run; - generated an orbital population plot; - reviewed the
existing Unity code and identified the gaps between the visualization
prototype and a scientific model.

## 1. Project Goal

This project aims to develop a scientifically grounded simulation of the
Kessler Syndrome: a cascading process in which collisions between
objects in Earth orbit generate debris and increase the likelihood of
subsequent collisions.

The intended system combines: - real orbital catalog data from TLE
sources; - SGP4 orbit propagation; - orbital mechanics calculations; -
close-approach and potential collision analysis; - collision and
fragmentation physics; - long-term population evolution; - interactive
3D visualization.

**Important:** this is currently a prototype. Successful execution and
agreement with simple analytical tests do not yet establish that the
complete model is a scientifically validated simulation.

## 2. Architecture: Python and Unity

### Python -- scientific computation backend

Python is being developed for: - TLE acquisition and storage; - TLE
syntax and checksum validation; - SGP4 state propagation; - osculating
orbital element calculations from state vectors; - orbital population
analysis; - future close-approach screening, TCA and collision
modelling.

### Unity -- existing visualization and dynamics prototype

Unity is not being discarded. The supplied Unity project already
contains: - a scene with Earth, camera and `SatelliteSpawner`; - an
Earth texture and satellite prefab; - loading of `satellites_state.json`
from `StreamingAssets`; - instantiation of catalogued objects; -
category-based object colouring; - a simplified Earth-centred gravity
model; - simulation time scaling and pause controls; - preliminary
spatial collision screening; - a test event triggered by the `C` key; -
an on-screen simulation statistics panel.

Python and Unity are **not yet integrated into one shared simulation
engine**. Unity currently loads its own JSON snapshot and propagates
object motion independently. The target architecture is for Python to
provide computed states and events, while Unity visualizes them.

## 3. Completed Python Work

### 3.1. Project structure and data pipeline

-   Organized raw data, processed data, metadata, source code, results
    and documentation directories.
-   Moved download and processing scripts into `src/data_pipeline/`.
-   Updated paths to resolve from the project root rather than the
    script directory.
-   Stored raw TLE groups under `data/raw/tle/`.

### 3.2. TLE acquisition and validation

Implemented: - downloading CelesTrak groups for active objects and
selected debris populations; - saving groups to separate TLE files; -
line-length and checksum validation; - NORAD identifier consistency
checks across TLE line pairs; - parsing through the `sgp4` package; -
propagation checks at each TLE epoch; - JSON and CSV validation reports.

Last reported validation run:

  Metric                               Result
  ---------------------------------- --------
  TLE files                                 5
  Records parsed                       18,582
  Valid records                        18,582
  Invalid records                           0
  Incomplete records                        0
  Warnings                                  6
  Duplicate NORAD IDs across files          0

The six warnings were deferred. These results apply to that particular
data snapshot and validation run; they do not guarantee that the TLE
data remain current.

### 3.3. SGP4 processing

`process_TLE_to_json.py` processes TLE records with SGP4 and writes
object positions and velocities to JSON. The state vectors are expressed
in the TEME reference frame, in kilometres and kilometres per second.

The latest reported orbital-mechanics run loaded 18,479 objects. The
difference from the 18,582 records in the validation report has not yet
been reconciled.

### 3.4. Orbital mechanics

The module under `src/orbital_mechanics/` implements vector operations
and calculations for: - specific angular momentum; - eccentricity vector
and eccentricity; - specific orbital energy; - semi-major axis; -
inclination; - RAAN; - argument of periapsis; - true anomaly; - perigee
and apogee altitude; - orbital period.

The elements are **instantaneous two-body osculating elements** derived
from SGP4 position and velocity vectors. They are not invariant orbital
parameters and should not be confused with mean TLE elements.

Circular-orbit test: - position: `(7000, 0, 0)` km; - velocity:
approximately `(0, 7.54605329, 0)` km/s; - semi-major axis:
approximately 7000 km; - eccentricity: approximately zero; -
inclination: 0 degrees; - perigee and apogee: approximately 629 km; -
orbital period: approximately 97.14 minutes.

An additional inclined elliptical test produced non-zero eccentricity
and inclination, but its selected initial conditions resulted in a
perigee below Earth's surface. It was a mathematical test, not a
physically viable satellite orbit.

Latest reported population processing:

  Metric                     Result
  ------------------------ --------
  Input objects              18,479
  Successfully processed     18,479
  Skipped                         0

Generated outputs: - `data/processed/orbital_elements.json` -
`results/figures/orbital_population_map.png`

The plot shows semi-major-axis altitude against inclination. It is not a
plot of instantaneous positions or spatial density.

## 4. Review of the Supplied Unity Project

See [docs/UNITY_AUDIT_EN.md](docs/UNITY_AUDIT_EN.md) for the detailed
code review.

Main conclusion: Unity is a useful visualization and simplified
orbital-dynamics prototype, but its current collision implementation is
not yet a physically reliable collision model.

Key limitations: - motion is independently propagated in Unity using
first-order semi-implicit Euler integration; - large time steps can
introduce significant numerical error; - collision screening checks
objects in the same spatial hash cell only and does not inspect
neighbouring cells; - detection checks instantaneous separation rather
than the minimum separation over the interval; - the 2 km threshold is a
close-approach/conjunction threshold, not physical contact between
objects; - collision handling relabels existing objects and applies
random velocity kicks, but does not create fragments or conserve
momentum and energy; - the `C` key triggers an artificial test event
rather than a computed orbital encounter; - the JSON bundled with Unity
contains 18,579 objects, different from the latest Python count of
18,479; the datasets need to be synchronized; - there is currently no
Python--Unity integration.

This does not make Unity useless. Its scene, Earth rendering, camera,
prefabs, population display, UI and event animation are a useful
foundation.

## 5. Current Scientific Limitations

1.  SGP4 outputs states in TEME. Reference frame, epoch and
    transformations must be made explicit before combining states or
    comparing them with other models.
2.  Position and velocity vectors used in one calculation must share a
    common epoch.
3.  Osculating elements derived from an instantaneous state must not be
    treated as identical to mean TLE elements.
4.  The two-body model does not automatically include J2, atmospheric
    drag, solar radiation pressure, or lunar and solar perturbations.
5.  Current object masses and radii assigned by category are estimates,
    not individual measurements.
6.  Collision probability, state uncertainty and a validated
    fragmentation model have not yet been implemented.

## 6. Development Roadmap

### Stage 4 -- Data consistency and population analysis

-   Reconcile record counts across validation output, Python JSON and
    Unity JSON.
-   Store data epoch and reference frame in metadata.
-   Inspect the population plot and produce altitude/inclination
    distributions.
-   Add unit tests for orbital element calculations.

### Stage 5 -- Collision Detection

-   Propagate candidate objects to a common epoch and reference frame.
-   Implement efficient candidate screening.
-   Calculate relative position and velocity.
-   Calculate Time of Closest Approach (TCA) and miss distance.
-   Account for motion between discrete steps to avoid missing fast
    encounters.
-   Distinguish close approaches from physical collisions.

### Stage 6 -- Collision Probability

-   Model position and velocity uncertainty.
-   Include object dimensions and encounter geometry.
-   Estimate collision probability only after defining an uncertainty
    model.

### Stage 7 -- Collision Physics and Fragmentation

-   Calculate relative impact velocity and collision energy.
-   Select a documented fragmentation model.
-   Generate fragment mass, size and velocity distributions.
-   Check mass and momentum budgets within the model assumptions.

### Stage 8 -- Debris Propagation

-   Propagate generated fragments.
-   Include appropriate perturbations and atmospheric drag.
-   Update the population over time.

### Stage 9 -- Kessler Cascade

-   Process subsequent collision events.
-   Track event history and fragment parentage.
-   Compare scenarios and perform sensitivity analysis.

### Stage 10 -- Unity Integration

-   Define one data exchange format.
-   Include stable NORAD identifiers and timestamps.
-   Unity should visualize backend-computed states and events rather
    than running an independent, divergent simulation.
-   Preserve Unity as the interactive 3D layer.

## 7. Immediate Next Objective

The next priorities are: 1. Reconcile and synchronize Python and Unity
input datasets. 2. Implement and test close-approach detection in
Python. 3. Validate the results against controlled cases. 4. Only then
feed computed conjunction and collision events into Unity.

**Progress criterion:** every module must have documented assumptions,
tests and known limitations. A successful animation or a clean program
run is not scientific validation.
