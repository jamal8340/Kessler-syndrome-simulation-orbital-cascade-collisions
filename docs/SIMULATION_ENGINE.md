# Simulation engine (`src/kessler_sim`)

The backend's physics engine. It has no networking; the Unity server (step C) is built on top of it.

## Running

```bash
uv sync                                            # installs the engine as a package
uv run pytest                                      # all tests (~5 s)
uv run python -m kessler_sim.demo scenarios        # list scenarios
uv run python -m kessler_sim.demo run iridium_cosmos_2009 --times 0,1h,1d,30d --plot
uv run python -m kessler_sim.demo run closest_approach_breakup
uv run python -m kessler_sim.demo approaches --duration 10m --top 10
```

`--plot` saves a figure to `results/figures/cloud_<scenario>.png` (or the given path).

## Server for Unity

```bash
uv run python -m kessler_sim.server --scenario catalogue          # ws://127.0.0.1:8765
uv run python -m kessler_sim.server --scenario iridium_cosmos_2009 --play --speed 600

# test client (another terminal), no Unity needed:
uv run python -m kessler_sim.client info          # session info and first state message
uv run python -m kessler_sim.client watch         # play 5 s, count state messages (~20/s)
uv run python -m kessler_sim.client demo          # closest approach -> collide -> seek across it
```

The message protocol is described in `../TRANSPORT_ARCHITECTURE.md` (in the parent `sims/` folder).

| Module | What it does |
|---|---|
| `server.py` | WebSocket server: streaming loop, commands, background jobs |
| `protocol.py` | Message encoding (JSON without infinity/NaN, binary state messages) |
| `earth.py` | Earth's rotation angle (Greenwich sidereal angle, as SGP4) |
| `client.py` | Command-line test client; also used by the tests |

## Modules

| Module | What it does |
|---|---|
| `orbits.py` | Orbit model: position/velocity of any object at any time t; state ↔ orbital elements |
| `state.py` | Object table: catalogue objects + fragments, IDs, parents, mass, size, orbit, lifespan |
| `clock.py` | Sim time, play / pause / speed / seek |
| `catalogue.py` | Loads the shared dataset; orbits from the raw TLEs (SGP4 mean elements at the epoch) |
| `breakup.py` | Collision → fragments (NASA breakup model from `collision_physics/fragmentation.py`) |
| `screening.py` | Close-approach search over a time window |
| `engine.py` | `Simulation`: the single entry point (state, breakups, close approaches, history) |
| `scenario.py`, `scenarios/*.json` | Scenario files: what to load, which breakups to schedule |
| `analysis.py` | Debris-cloud statistics over time |
| `demo.py` | Command-line demo |

## Using it from Python

```python
from kessler_sim import load_scenario

sim = load_scenario("catalogue")                  # 18,558 objects, t = 0 is the dataset epoch
snap = sim.state(3600)                            # ids, positions [km], velocities [km/s] at t = 1 h
found = sim.find_close_approaches(0, 600)         # encounters in the first 10 minutes
event = sim.collide_approach(found[0])            # force the closest pair to collide at its TCA
print(event.summary())                            # fragment counts, masses, energy, ...
sim.state(event.t - 1)                            # before the breakup: parents are back
```

## Model and simplifications

- **Orbits from a formula.** Two-body motion plus the secular (long-term) J2 drift, i.e. the slow
  turning of the orbit caused by Earth's flattening. No time-stepping: any time can be computed
  directly, so seeking is free and runs are reproducible.
- **Catalogue orbits** are SGP4's mean elements at the dataset epoch, plus a per-object decay rate
  fitted from SGP4 over one day (the only drag effect).
- **Fragment orbits** come from their position and velocity at the breakup.
- **Lifespans.** Each object exists in `[t_created, t_destroyed)`. A breakup destroys both parents
  and creates fragments at the same instant. Seeking to before a breakup shows the parents again.
- **Rewriting history.** Triggering a breakup at time t discards all later breakups. Fragment IDs
  are deterministic (`1,000,000 + event × 100,000 + k`), so replaying the same breakup gives the
  same fragments.
- **Two debris clouds.** In a catastrophic collision each fragment is assigned to one parent and
  flies with that parent's velocity plus its kick. Each parent's momentum is conserved exactly, and
  the debris keeps the parent's orbital plane, as observed for real collisions.
- **Non-catastrophic collision.** The larger object survives as a `[REMNANT]` object.
- **Collision point.** All fragments start at the heavier object's position. For pairs from the
  close-approach list this is exact to a few km. For arbitrary "what if" pairs, the event's
  `separation_km` shows how far apart the objects really were.
- **Re-entry.** Fragments with perigee below 150 km (configurable) or on escape trajectories are
  counted but not added.
- **Fragment cap.** At most 5,000 fragments per breakup (configurable); above that, the minimum
  fragment size is raised instead.

## Accuracy (vs SGP4, low Earth orbit)

| Time after epoch | Median error | 95th percentile |
|---|---|---|
| 0 | ~10 km | ~17 km |
| 1 day | ~10 km | ~17 km |
| 3 days | ~10 km | ~20 km |
| 7 days | ~11 km | ~100 km |

The ~10 km comes from the short-period J2 terms we ignore. Accuracy is much worse (up to
thousands of km) for 15 very elongated deep-space orbits (e.g. the MMS satellites), where SGP4
adds Moon and Sun effects; they are outside the project's scope.

## Known limitations

- About 5% of fragments get kicks above ~1 km/s from the tail of the breakup model's velocity
  distribution and end up on very elongated orbits (apogee above 2,000 km).
- Masses and sizes of catalogue objects are per-category estimates.
- No drag on fragments other than the re-entry cut-off; no Moon/Sun perturbations.
- See also the known limitation recorded by the `xfail` test in `tests/test_fragmentation.py`.
