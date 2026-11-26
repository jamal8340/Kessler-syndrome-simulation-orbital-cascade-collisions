# Kessler Syndrome Simulation – backend

Simulation of collisions between objects in Earth orbit: how two objects break into fragments, what is left of them, where the pieces fly and how the debris cloud spreads over time. It uses the real catalogue of ~18,500 tracked satellites and debris pieces (CelesTrak data) and NASA's standard breakup model.

This repository is the **backend**: all physics runs here. The 3D visualisation is the Unity project `KesslerSimulation` (see its `README.md`), which connects to this backend's server. Without Unity, the backend can be used from the command line.

## What it can do

- Show where every catalogued object is at any moment, forwards or backwards in time.
- Find pairs of objects that will pass close to each other.
- Collide any two objects; the resulting fragments become new objects on their own orbits (two debris clouds, one per parent); fragments that would re-enter immediately are removed.
- Rewind and replay: before a collision the original objects are back; replaying gives the same fragments.
- Replay two real events: the Iridium 33 × Cosmos 2251 collision (2009) and the Fengyun‑1C anti-satellite test (2007), and compare the simulated debris with the real debris.
- Serve all of this live to the Unity visualisation.

## Requirements and setup

- [uv](https://docs.astral.sh/uv/) – it installs the right Python version (3.14) and all dependencies automatically on the first `uv run`.
- The data needed to run is included in the repository (`data/`).

```bash
uv sync          # optional: installs everything up front
```

All commands below are run from this repository's folder.

## Long-term drag (F1, experimental)

The optional orbit-averaged drag model preserves orbital history and re-entry times. It is not enabled in the existing live server; it is intended for the future cascade batch runner. Run the two reference fragments or the full catalogue without collisions:

```bash
uv run python -m kessler_sim.drag_demo --years 50
uv run python -m kessler_sim.drag_demo --scenario catalogue --years 10
uv run python -m kessler_sim.drag_demo --scenario catalogue --years 10 --no-stationkeeping
uv run pytest tests/test_drag.py tests/test_stationkeeping.py -q
```

By default, intact active satellites ideally compensate drag (constant orbital size and eccentricity); fragments, debris and remnants remain passive. Use `--no-stationkeeping` for the original all-passive comparison. This assumes indefinite operation, not realistic mission lifetimes.

See `docs/ATMOSPHERIC_DRAG.md` for assumptions, the Python API, and limitations.

## Quick start with Unity

```bash
uv run python -m kessler_sim.server --play --speed 600
```

Then open the Unity project and press Play (see the Unity `README.md`). The server prints `listening on ws://127.0.0.1:8765` when it is ready; it keeps one simulation for everyone who connects. Stop it with Ctrl+C; to switch scenario or start over, restart it.

### Server options

| Option | Default | Meaning |
|---|---|---|
| `--scenario NAME` | `catalogue` | Which scenario to load (see below) |
| `--play` | off | Start with the clock running |
| `--speed S` | 60 | Simulated seconds per real second |
| `--port P` | 8765 | Port Unity connects to |
| `--host H` | 127.0.0.1 | Network address to listen on |
| `--rate R` | 20 | Position updates per second sent to Unity |

## Scenarios

| Scenario | What it shows |
|---|---|
| `catalogue` | All 18,558 catalogued objects at the dataset's moment (2026‑11‑12 22:50 UTC); no collisions – a playground for your own collisions |
| `closest_approach_breakup` | Full catalogue; the closest pair passing within the first 10 minutes collides (currently STARLINK‑6135 × SIGI, ~766 fragments) |
| `iridium_cosmos_2009` | Replay of the 2009 collision of Iridium 33 and Cosmos 2251 (only these two satellites and their debris); collision at t = 600 s |
| `fengyun_1c_2007` | Replay of the 2007 anti-satellite test: an interceptor rises from the ground and destroys the Fengyun‑1C weather satellite; collision at t = 600 s |

`uv run python -m kessler_sim.demo scenarios` lists them. The historical scenarios use idealised orbits with the real inclinations, altitudes and collision geometry; see `docs/VALIDATION.md`.

## Using it without Unity

**Demo** – what happened in a scenario and how the debris cloud spreads:
```bash
uv run python -m kessler_sim.demo run iridium_cosmos_2009 --times 0,10m,1h,1d,30d --plot
uv run python -m kessler_sim.demo approaches --duration 10m --top 10     # close passes in the catalogue
```
`--plot` saves a figure to `results/figures/cloud_<scenario>.png`.

**Test client** – checks a running server without Unity (in a second terminal):
```bash
uv run python -m kessler_sim.client info     # what the server sends on connect
uv run python -m kessler_sim.client watch    # play for 5 s and count the updates
uv run python -m kessler_sim.client demo     # closest pass -> collide -> jump before and after
```

**Validation** – compares the simulated debris of the two historical events with the real debris still in orbit today:
```bash
uv run python -m kessler_sim.validation                              # both events
uv run python -m kessler_sim.validation --scenario fengyun_1c_2007   # one of them
```
It writes figures to `results/figures/validation_*.png` and tables to `results/validation_*.json` / `.md`. The results and their interpretation are in `docs/VALIDATION.md`.

## Updating the data

The catalogue is a snapshot downloaded from CelesTrak (currently 2026‑11‑12). To use newer data, run these steps in order:

| Step | Command | Result |
|---|---|---|
| 1. Download | `uv run python src/data_pipeline/TLE_data_download.py` | Active satellites and four debris groups in `data/raw/tle/` (needs internet) |
| 2. Check | `uv run python src/data_pipeline/validate_tle.py` | Report in `data/metadata/tle_validation_report.json` / `_details.csv` |
| 3. Process | `uv run python src/data_pipeline/process_TLE_to_json.py --epoch 2026-11-12T22:50:04+00:00` | The shared dataset `data/processed/satellites_state.json` and its summary; the dataset's moment (epoch) is the given time, or now if `--epoch` is left out |
| 4. Unity snapshot | `uv run python src/data_pipeline/export_to_unity.py` | Copies the dataset into the Unity project for its offline mode (if the previous file differs, it is backed up in `KesslerSimulation/DatasetBackups/`) |

Afterwards restart the server. Rerunning step 3 with the same `--epoch` reproduces the same dataset exactly.

## Tests

```bash
uv run pytest        # ~145 tests, about 12 s
```
One test is marked as an expected failure: it documents a known limitation of the breakup model (rare weak impacts).

## Accuracy and limitations

- Positions agree with the SGP4 model to about 10 km in low Earth orbit over several days (much worse for a few very elongated deep-space orbits, which are out of scope).
- Collisions only happen when triggered (by you, by a scenario, or from the close-approach list); nothing collides by itself yet. A long-term cascade with statistically sampled collisions at real sizes is planned.
- Masses and sizes of catalogue objects are estimates by object type.
- Compared with the real 2009 and 2007 debris: orbital planes match, but fragment kicks are too strong (clouds spread wider than in reality), fragments have no atmospheric drag, and the Fengyun‑1C test produces too few fragments. Details in `docs/VALIDATION.md`.

## Older tools

These scripts from earlier project stages still work; the first two are superseded by the `kessler_sim` package above:

| Script | Status |
|---|---|
| `src/collision_detection/close_approach.py` (10-minute close-approach scan, writes `data/processed/close_approaches.json`) | superseded by `kessler_sim.demo approaches` and the Unity close-pass panel |
| `src/collision_physics/run_breakups.py` (breakup statistics for close approaches) | superseded by collisions in `kessler_sim` (demo, server, Unity) |
| `src/orbital_mechanics/orbital_elements.py` (orbital elements and the altitude/inclination plot `results/figures/orbital_population_map.png`) | not superseded – still the only overview plot of the whole population |

## Further documentation

| Document | Contents |
|---|---|
| `docs/SIMULATION_ENGINE.md` | How the simulation engine works, its model and accuracy |
| `docs/VALIDATION.md` | Comparison with the real debris of the 2009 and 2007 events |
| `../KesslerSimulation/README.md` | The Unity visualisation: controls, panels, watching collisions |
| `docs/UNITY_AUDIT_EN.md.old`, `docs/UNITY_AUDIT_PL.md.old` | Outdated: review of the old Unity prototype, before it was connected to this backend |

A Polish version of this file is `README_PL.md`.
