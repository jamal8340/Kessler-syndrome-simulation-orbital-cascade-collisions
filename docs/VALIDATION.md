# Validation: Iridium 33 × Cosmos 2251

```bash
uv run python -m kessler_sim.validation
```

Compares the simulated debris of the `iridium_cosmos_2009` scenario with the real debris of that collision, using only local data. The simulation is **not** adjusted to the results; differences are recorded as accepted limitations.

| What | Simulated | Real |
|---|---|---|
| Source | `iridium_cosmos_2009` scenario | CelesTrak TLE groups in `data/raw/tle/` (`iridium-33-debris`, `cosmos-2251-debris`) |
| Moment | right after the collision (2009‑02‑10 16:56 UTC) | fragments still in orbit at the TLE epochs (late 2026, ~17 years later) |

Outputs: `results/figures/validation_inclination.png`, `validation_gabbard.png`, `validation_perigee_apogee.png`, and the numbers in `results/validation_iridium_cosmos.json` / `.md`.

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

## Tests

`tests/test_validation.py` checks the data loading, that the mean inclinations agree within 0.2°, and records the accepted limitation (simulated spread wider than real) so that a change in the model's behaviour is noticed.
