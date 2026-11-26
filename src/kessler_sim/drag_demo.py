"""F1 smoke test and drag-only benchmark; no statistical collisions or Unity needed."""
from __future__ import annotations

import argparse
import math
import time

import numpy as np

from . import DragConfig, Simulation, load_scenario, orbits

YEAR = 365.25 * 86400


def reference_fragments():
    sim = Simulation(name="drag_reference")
    a = orbits.R_MEAN + np.array([400., 800.])
    z = np.zeros(2)
    sim.table.add([1, 2], ["400 km", "800 km"], "fragment", 1., .1,
                  orbits.Elements(a, z, z + .7, z, z, z), 0., am_m2_per_kg=.02)
    return sim


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", default="reference", help="reference (two fragments) or existing scenario name")
    parser.add_argument("--years", type=float, default=50.)
    parser.add_argument("--density-scale", type=float, default=1.)
    parser.add_argument("--no-stationkeeping", action="store_true",
                        help="treat all objects as passive (comparison run)")
    args = parser.parse_args()
    if not math.isfinite(args.years) or args.years <= 0:
        parser.error("--years must be finite and positive")
    sim = reference_fragments() if args.scenario == "reference" else load_scenario(args.scenario)
    start = sim.events[-1].t if sim.events else 0.
    before = sim.state(start).r_km.copy()
    begin = time.perf_counter()
    sim.enable_drag(DragConfig(density_scale=args.density_scale,
                               stationkeeping=not args.no_stationkeeping,
                               reentry_altitude_km=sim.config.reentry_altitude_km), t=start)
    end = start + args.years * YEAR
    sim.advance_drag(end)
    elapsed = time.perf_counter() - begin
    print(f"Drag-only run: {args.scenario}, {args.years:g} years, {len(sim.table)} historical objects")
    print(f"Station keeping: {sim.drag_config.stationkeeping}")
    print(f"Alive: {len(sim.state(end).ids)}, drag re-entries: {len(sim.reentries)}")
    alive = sim.table.alive_mask(end)
    for category in np.unique(sim.table.category):
        print(f"  {category}: {np.count_nonzero(alive & (sim.table.category == category))}")
    print(f"Segments: {sim.table.orbit_history.segment_count}; allocated history buffers: "
          f"{sim.table.orbit_history.nbytes / 2**20:.2f} MiB; wall time: {elapsed:.3f} s")
    if args.scenario == "reference":
        for oid in (1, 2):
            info = sim.object_info(oid, end)
            reentry = sim.reentries.get(oid)
            if reentry is None:
                print(f"{info['name']}: alive, perigee {info['orbit']['perigee_alt_km']:.2f} km")
            else:
                print(f"{info['name']}: re-entry after {(reentry-start)/86400:.2f} days")
        np.testing.assert_allclose(sim.state(start).r_km, before, atol=1e-8)
        print("Rewind check: original positions restored")


if __name__ == "__main__":
    main()
