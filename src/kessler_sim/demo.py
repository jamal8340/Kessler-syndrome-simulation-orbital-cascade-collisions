"""Command-line demo of the simulation engine (no Unity needed).

    uv run python -m kessler_sim.demo scenarios
    uv run python -m kessler_sim.demo run iridium_cosmos_2009 --times 0,1h,6h,1d,3d --plot
    uv run python -m kessler_sim.demo run closest_approach_breakup
    uv run python -m kessler_sim.demo approaches --duration 600 --top 10
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

from . import orbits
from .analysis import cloud_stats
from .scenario import available_scenarios, load_scenario

FIGURES_DIR = Path(__file__).resolve().parents[2] / "results" / "figures"
_UNITS = {"s": 1, "m": 60, "h": 3600, "d": 86400}


def parse_duration(text: str) -> float:
    m = re.fullmatch(r"\s*([0-9.]+)\s*([smhd]?)\s*", text)
    if not m:
        raise argparse.ArgumentTypeError(f"bad duration '{text}' (use e.g. 90, 10m, 6h, 3d)")
    return float(m.group(1)) * _UNITS[m.group(2) or "s"]


def fmt_dt(seconds: float) -> str:
    for unit, size in (("d", 86400), ("h", 3600), ("m", 60)):
        if abs(seconds) >= size:
            return f"{seconds / size:g}{unit}"
    return f"{seconds:g}s"


def print_event(ev):
    s = ev.summary()
    print(f"\nBreakup #{s['event_id']} at t = {s['t']:.1f} s: {s['name_a']} x {s['name_b']}")
    print(f"  impact speed {s['v_rel_kms']:.2f} km/s at {s['altitude_km']:.0f} km altitude "
          f"(objects {s['separation_km']:.2f} km apart)")
    print(f"  {'catastrophic' if s['catastrophic'] else 'non-catastrophic'} (EMR {s['emr_jg']:.0f} J/g), "
          f"masses {s['mass_a_kg']:.0f} kg + {s['mass_b_kg']:.0f} kg")
    print(f"  fragments >= {s['l_min_m']:.2f} m: {s['n_fragments']} in orbit "
          f"(model {s['n_expected']}; {s['reentered_count']} re-entered, {s['escaped_count']} escaped)")
    print(f"  per parent: {s['fragments_per_parent']}")
    print(f"  mass: {s['tracked_mass_kg']:.0f} kg tracked, {s['untracked_mass_kg']:.0f} kg untracked, "
          f"{s['lost_mass_kg']:.0f} kg re-entered/escaped")
    if s["survivor_id"] is not None:
        print(f"  survivor: object {s['survivor_id']}")
    if s["largest_fragment_mass_kg"] is not None:
        print(f"  largest fragment: {s['largest_fragment_mass_kg']:.1f} kg, {s['largest_fragment_lc_m']:.2f} m")


def print_cloud_table(sim, ev, offsets):
    print(f"\n  Cloud of breakup #{ev.event_id} over time:")
    print(f"  {'time after':>10} {'alive':>6} {'alt p5':>8} {'alt p95':>8} {'arc p95':>8}  "
          f"per parent: inclination mean±std, RAAN std")
    for dt in offsets:
        c = cloud_stats(sim, ev.event_id, ev.t + dt)
        if not c["n_alive"]:
            print(f"  {fmt_dt(dt):>10} {0:>6}")
            continue
        per = "; ".join(f"{sim.table.names[sim.table.row(pid)]}: {p['i_mean_deg']:.1f}±{p['i_std_deg']:.1f}°, "
                        f"{p['raan_std_deg']:.2f}°" for pid, p in c["per_parent"].items())
        print(f"  {fmt_dt(dt):>10} {c['n_alive']:>6} {c['altitude_p5_km']:>6.0f}km {c['altitude_p95_km']:>6.0f}km "
              f"{c['arc_spread_p95_deg']:>7.1f}°  {per}")


def plot_clouds(sim, offsets, path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tb = sim.table
    fig, axes = plt.subplots(len(sim.events), 2, figsize=(13, 5.5 * len(sim.events)), squeeze=False)
    for ax_row, ev in zip(axes, sim.events):
        ax_pos, ax_gab = ax_row
        rows = sim.fragment_rows(ev.event_id)
        origin = {int(i): int(o) for i, o in zip(ev.fragment_ids, ev.fragment_origin_ids)}
        colours = {ev.id_a: "#1f77b4", ev.id_b: "#d62728"}
        theta = np.linspace(0, 2 * np.pi, 200)
        ax_pos.plot(orbits.R_MEAN * np.cos(theta), orbits.R_MEAN * np.sin(theta), color="0.6", lw=1)
        for dt, alpha in zip(offsets, np.linspace(0.25, 1.0, len(offsets))):
            alive = rows[tb.alive_mask(ev.t + dt)[rows]]
            r, _ = tb.states(ev.t + dt, alive)
            c = [colours.get(origin.get(int(i)), "k") for i in tb.ids[alive]]
            ax_pos.scatter(r[:, 0], r[:, 1], s=2, c=c, alpha=alpha)
        ax_pos.set_aspect("equal")
        lim = 1.6 * orbits.R_MEAN
        ax_pos.set_xlim(-lim, lim)
        ax_pos.set_ylim(-lim, lim)
        ax_pos.set_title(f"Breakup #{ev.event_id}: fragments at t+{', '.join(fmt_dt(d) for d in offsets)} (x-y plane)")
        ax_pos.set_xlabel("x [km]")
        ax_pos.set_ylabel("y [km]")

        el = tb.elements_at(ev.t, rows)
        per = orbits.period(el) / 60.0
        c = [colours.get(origin.get(int(i)), "k") for i in tb.ids[rows]]
        ax_gab.scatter(per, orbits.apogee_radius(el) - orbits.R_MEAN, s=3, c=c, marker="^", label="apogee")
        ax_gab.scatter(per, orbits.perigee_radius(el) - orbits.R_MEAN, s=3, c=c, marker="v", label="perigee")
        ax_gab.set_title(f"Gabbard diagram (blue: {ev.name_a}, red: {ev.name_b})")
        ax_gab.set_xlim(np.percentile(per, 1) - 1, np.percentile(per, 97) + 1)
        ax_gab.set_ylim(max(0, np.percentile(orbits.perigee_radius(el) - orbits.R_MEAN, 1) - 50),
                        np.percentile(orbits.apogee_radius(el) - orbits.R_MEAN, 97) + 50)
        ax_gab.set_xlabel("orbital period [min]")
        ax_gab.set_ylabel("altitude [km]")
        ax_gab.grid(alpha=0.3)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"\nPlot saved: {path}")


def cmd_run(args):
    sim = load_scenario(args.scenario)
    print(f"Scenario '{sim.name}': {getattr(sim, 'description', '')}")
    print(f"Epoch {sim.epoch_utc:%Y-%m-%d %H:%M:%S} UTC, {len(sim.table)} objects in the table, "
          f"{len(sim.events)} breakup(s)")
    offsets = [parse_duration(x) for x in args.times.split(",")]
    for ev in sim.events:
        print_event(ev)
        print_cloud_table(sim, ev, offsets)
        high = cloud_stats(sim, ev.event_id, ev.t)["high_orbit_count"]
        print(f"  ({high} fragments have apogees above 2000 km: tail of the model's velocity-kick distribution)")
    if args.plot and sim.events:
        plot_clouds(sim, offsets, Path(args.plot) if isinstance(args.plot, str) else
                    FIGURES_DIR / f"cloud_{sim.name}.png")


def cmd_approaches(args):
    sim = load_scenario(args.scenario)
    found = sim.find_close_approaches(args.start, args.duration)
    print(f"{len(found)} close approaches (<= 10 km) between t = {args.start:g} s and "
          f"t = {args.start + args.duration:g} s:")
    for e in found[:args.top]:
        print(f"  {e['name_1'][:24]:24} x {e['name_2'][:24]:24} {e['miss_distance_km']:7.3f} km  "
              f"TCA {e['tca_s']:7.1f} s  {e['relative_speed_kms']:5.2f} km/s  {e['event_type']}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Kessler simulation engine demo")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("scenarios", help="list available scenarios")
    run = sub.add_parser("run", help="load a scenario and report its breakups and debris clouds")
    run.add_argument("scenario")
    run.add_argument("--times", default="0,10m,1h,6h,1d,3d", help="time offsets after each breakup")
    run.add_argument("--plot", nargs="?", const=True, default=False,
                     help="save a figure (optionally to the given path)")
    appr = sub.add_parser("approaches", help="list close approaches in a time window")
    appr.add_argument("--scenario", default="catalogue")
    appr.add_argument("--start", type=parse_duration, default=0.0)
    appr.add_argument("--duration", type=parse_duration, default=600.0)
    appr.add_argument("--top", type=int, default=10)
    args = ap.parse_args(argv)

    if args.cmd == "scenarios":
        for name in available_scenarios():
            print(name)
    elif args.cmd == "run":
        cmd_run(args)
    else:
        cmd_approaches(args)


if __name__ == "__main__":
    main()
