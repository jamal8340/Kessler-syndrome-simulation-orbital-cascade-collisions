# Validation against real debris

```bash
uv run python -m kessler_sim.validation                              # both scenarios
uv run python -m kessler_sim.validation --scenario fengyun_1c_2007   # one of them
```

Compares simulated debris with the real debris of the same event, using only local data, for two historical scenarios:
1. **Iridium 33 × Cosmos 2251** (`iridium_cosmos_2009`), the 2009 accidental collision;
2. **Fengyun‑1C anti-satellite test** (`fengyun_1c_2007`), see the section further below.

The simulation is **not** adjusted to the results; differences are recorded as accepted limitations.

# 1. Iridium 33 × Cosmos 2251

| What | Simulated | Real |
|---|---|---|
| Source | `iridium_cosmos_2009` scenario | CelesTrak TLE groups in `data/raw/tle/` (`iridium-33-debris`, `cosmos-2251-debris`) |
| Moment | right after the collision (2009‑02‑10 16:56 UTC) | fragments still in orbit at the TLE epochs (late 2026, ~17 years later) |

Outputs: `results/figures/validation_iridium_cosmos_inclination.png`, `…_gabbard.png`, `…_perigee_apogee.png`, and the numbers in `results/validation_iridium_cosmos.json` / `.md`.

## Results

| Quantity | Iridium 33: simulated | Iridium 33: real | Cosmos 2251: simulated | Cosmos 2251: real |
|---|---|---|---|---|
| Fragments | 536 | 109 | 561 | 581 |
| Inclination mean | 86.38° | 86.33° | 74.04° | 74.01° |
| Inclination std | 0.34° | 0.08° | 0.15° | 0.07° |
| Perigee 5–95 % | 325–789 km | 527–760 km | 303–789 km | 514–769 km |
| Apogee 5–95 % | 789–2221 km | 541–922 km | 789–2012 km | 584–1324 km |
| Apogee above 1500 km | 10 % | 0 % | 10 % | 1 % |

## Interpretation

**Confirmed – orbital planes.** Inclination hardly changes over the years, so it is the most reliable comparison. The mean inclinations agree within 0.05°: fragments keep their parent's orbital plane, which confirms the two-cloud breakup design (each fragment flies with its parent's velocity plus its kick).

**Explained by missing drag – perigees, apogees below ~800 km, counts.** The simulation has no atmospheric drag for fragments (only an immediate re-entry cut at 150 km). In the 17 years before the real data was recorded:
- fragments with low perigees decayed, so the real debris has no perigees below ~450 km while the simulated cloud reaches ~300 km;
- drag lowered and circularised the remaining orbits, which is why the real Gabbard diagram's "X" is shifted below the collision altitude (789 km);
- most fragments re-entered, especially Iridium's (109 real vs 536 simulated); the closeness of the Cosmos counts (581 vs 561) is a coincidence.

**Accepted limitation – kicks too strong.** The simulated cloud is wider than the real one in two ways that drag does not explain:
- the inclination spread is 2–4× larger (0.34° vs 0.08°, 0.15° vs 0.07°);
- 10 % of simulated fragments have apogees above 1500 km, against 0–1 % of real ones. Fragments with perigees near 789 km and high apogees decay slowly, so drag cannot have removed most of them.

Both point to fragment velocity kicks that are too large: the long tail of the NASA breakup model's kick distribution as implemented here (median 87 m/s, but 5 % above ~530 m/s and a maximum near 3 km/s). This is accepted as a limitation of the model and not corrected.

**Not compared.** Fragment sizes and masses (the TLEs do not contain them), and the exact 2009 orbits (the scenario uses idealised circular orbits with the real inclinations, altitude and impact geometry).

# 2. Fengyun‑1C anti-satellite test (2007)

| What | Simulated | Real |
|---|---|---|
| Source | `fengyun_1c_2007` scenario | CelesTrak TLE group `fengyun-1c-debris` in `data/raw/tle/` |
| Moment | right after the impact (2007‑01‑11, ~22:26 UTC) | fragments still in orbit at the TLE epochs (late 2026, ~20 years later) |

**The scenario.** A kinetic kill vehicle on a suborbital trajectory hits the Fengyun‑1C weather satellite (sun-synchronous, ~865 km, 98.6°) at ~8 km/s. Fengyun‑1C is on an idealised circular orbit with its real inclination and altitude; the interceptor is defined by its position and velocity at impact (steep climb) and exists only from its launch, 277 s before impact. Masses are approximate (Fengyun‑1C ~880 kg; interceptor ~600 kg, a public estimate); the time of day is approximate.

Only Fengyun‑1C's own debris is compared: the real data contains no interceptor debris, and in the simulation all interceptor fragments inherit its suborbital path and fall back immediately (661 fragments re-entered).

Outputs: `results/figures/validation_fengyun_1c_inclination.png`, `…_gabbard.png`, `…_perigee_apogee.png`, and `results/validation_fengyun_1c.json` / `.md`.

## Results

| Quantity | Fengyun‑1C: simulated | Fengyun‑1C: real |
|---|---|---|
| Fragments | 562 | 1987 |
| Inclination mean | 98.68° | 98.92° |
| Inclination std | 1.62° | 0.73° |
| Perigee 5–95 % | 408–865 km | 554–851 km |
| Apogee 5–95 % | 865–2373 km | 686–1300 km |
| Apogee above 1500 km | 13 % | 2 % |

## Interpretation

**Orbital plane – roughly confirmed.** The mean inclinations differ by 0.24°, more than for Iridium × Cosmos (0.05°). Possible causes, not verified here: the scenario's input inclination (98.6°) is approximate, and sun-synchronous orbits slowly change inclination over 20 years under the Sun's gravity, which the model does not include.

**Accepted limitation – kicks too strong (same as above).** The inclination spread is 2.2× wider than real (1.62° vs 0.73°), and 13 % of simulated apogees are above 1500 km against 2 % real. The real cloud does have a long high-apogee tail too (up to ~2,600 km, visible in the Gabbard diagram), stronger than for Iridium × Cosmos, but the simulated tail is heavier.

**Accepted limitation – too few fragments.** The model leaves 562 fragments (≥ 10 cm) in orbit, while 1,987 are still in orbit in reality after ~20 years of drag. So the model under-predicts this event by more than 3×, even against the survivors. Two parts of the model contribute: the fragment count formula of the NASA breakup model for the given masses (1,224 fragments in total), and the split of fragments between the two objects by their mass budgets, which here gives the interceptor 610 of the 1,224 fragments, all of which fall back (of Fengyun‑1C's 614 fragments, 562 stay in orbit).

**Explained by missing drag.** As for Iridium × Cosmos: no low perigees in the real data, the real Gabbard "X" sits lower, and the count comparison is only a lower bound.

## Tests

`tests/test_validation.py` checks the data loading of both scenarios, that the mean inclinations agree (within 0.2° for Iridium × Cosmos, 0.5° for Fengyun‑1C), and records the accepted limitations (simulated spread wider than real; interceptor fragments all fall back) so that a change in the model's behaviour is noticed.
