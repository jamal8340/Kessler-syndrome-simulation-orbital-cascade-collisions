# F1 — long-term atmospheric drag

## Model and scope

The optional drag model is intended for the statistical cascade (F2 onward), not precise orbit prediction. Existing scenarios and the server retain the legacy model by default. Enable drag only after the last preloaded breakup. No collision radii, fragment model, or protocol format are changed.

For a spherical, non-rotating atmosphere, define B = Cd A/m and drag acceleration = -0.5 rho B |v| v. Cd = 2.2 by default. Fragment A/m comes from the breakup model; when absent, use pi r^2/m from the existing category estimates. For a bound Kepler orbit with semimajor axis a, eccentricity e, and true anomaly nu, energy and angular momentum give:

- da/dt = -B rho v^3 a^2 / mu (all quantities in SI, converted back to km/s in code).
- de/dt = -B rho v (e + cos nu).

The secular rates are time-averaged using 48-point Gauss-Legendre quadrature over eccentric anomaly E in [0, pi], exploiting symmetry and weighting by dM/dE = 1 - e cos E. Circular orbits stay circular. A constant nominal density profile is log-linearly interpolated between tabulated heights from 150–1000 km, clamped below 150 km, extrapolated above 1000 km with the final interval slope. The tabulated densities follow Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed. (2013); independently accessible numerical reference: https://github.com/lcx366/ATMOS/blob/master/pyatmos/standardatmos/expo.py . The interpolation is our own continuous approximation. The `density_scale` parameter allows sensitivity analysis, not fitting to historical results.

Only bound orbits with perigee <= 2000 km are subject to the new model (including eccentric LEO-crossing orbits). Orbits already at or below the re-entry threshold are marked destroyed at activation/birth. Above the scope, passive objects retain the original orbital model. Intact active objects are handled separately by station keeping. Managed objects replace the legacy catalogue `ndot` at activation, preserving their initial phase. The initial elements and t_ref are never overwritten.

## Ideal station keeping (default)

`DragConfig(stationkeeping=True)` assumes that intact objects with category `active` and `parent_event == -1` compensate all atmospheric drag indefinitely. Their semimajor axis and eccentricity stay constant, legacy `ndot` is disabled at activation, and J2 phase propagation continues. A single open-ended orbital segment preserves the earlier history without repeated computation. The re-entry threshold still takes precedence. This policy also covers active objects outside the drag altitude range; passive objects outside that range retain legacy propagation.

Fragments and remnants created by collisions are always passive, even if a remnant retains the `active` display category. `DragConfig(stationkeeping=False)` restores the original all-passive F1 experiment. No actual manoeuvres, propulsion inventory, failure rates, retirements or replacement launches are modelled. Category `active` is a proxy for operational capability, not verified propulsion status. This is an idealized maintained-population experiment, not a 50-year forecast.

## Segments and API

`DragConfig` defaults: 90-day maximum segment, 5 km maximum semimajor axis change, maximum eccentricity change 0.002, and 150 km re-entry threshold. Each packed segment stores start/end time, six mean elements, two constant decay rates, and three secular J2 phase rates at the segment midpoint. Midpoint drag rates are estimated from the preceding rates; the step is shortened as needed. Elements vary continuously across boundaries. The perigee crossing is located by bisection within the linear segment. Position is continuous; velocity is the instantaneous Kepler velocity, not a derivative of the secular interpolation. This is not a high-fidelity trajectory integrator.

```python
from kessler_sim import Simulation, DragConfig, BreakupParams

sim = Simulation.from_catalogue()
sim.enable_drag(DragConfig(), t=0.)
sim.advance_drag(3600.)  # before collision sampling at 3600 s
# event = sim.trigger_breakup(id_a, id_b, 3600., BreakupParams())
sim.advance_drag(86400.)
old = sim.state(0.)      # read-only seek
sim.truncate_after(3600.) # discard future history and re-entries
```

`advance_drag(t)` advances the committed horizon monotonically and returns newly committed `(object_id, time)` re-entries in chronological order. `sim.reentries` maps IDs to committed times (not immediate losses in breakup events). New fragments get segments on the next advance; segments can predict beyond the committed horizon, but predictions do not destroy objects. Always advance drag before using a later time for physics; querying beyond the committed horizon is only a preview and may extrapolate the final planned orbit. `trigger_breakup` rejects future times beyond that horizon. Rewinding with `state`/`object_info` is read-only, whereas `truncate_after` removes later segments and re-entries, revives objects, and removes later breakups. Rewriting before activation disables drag; enable it again explicitly. Drag segment storage is indexed by stable object ID, with vectorized latest-segment access and compact historical arrays. `ObjectTable.elements_at(t, rows)` is the time-aware element accessor.

The server is not yet wired to advance drag or constrain its clock to a computed horizon. Replay/persistence and the batch runner are F3, Unity controls F5. Do not enable drag on the existing live server until those stages are integrated.

## Verification and limitations

```bash
uv run pytest tests/test_drag.py tests/test_stationkeeping.py -q
uv run pytest
uv run python -m kessler_sim.drag_demo --years 50
uv run python -m kessler_sim.drag_demo --scenario catalogue --years 10
uv run python -m kessler_sim.drag_demo --scenario catalogue --years 10 --no-stationkeeping
```

Reference fragments: circular at 400 and 800 km, A/m=0.02 m²/kg, Cd=2.2. Nominal 400 km re-entry around 70 days; 800 km survives 50 years. Actual lifetimes depend strongly on solar/geomagnetic activity, altitude, attitude and A/m. This static atmosphere neglects those variations, winds, atmospheric rotation, actual station-keeping manoeuvres, mission lifetimes, launches, and solar radiation pressure. The catalogue estimates can cause substantial attrition when no station keeping or replacements are modelled; do not interpret drag-only population decline as a forecast. The quadrature and segment tolerances should be checked for sensitivity, especially very eccentric orbits and rapid decay near re-entry where orbit averaging is least reliable. Long-term orbital phase is approximate. No empirical drag validation or cascade collision-rate validation is claimed.
