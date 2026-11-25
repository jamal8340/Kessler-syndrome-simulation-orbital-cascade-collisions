"""Step E: compare simulated debris with the real debris of the same event.

    uv run python -m kessler_sim.validation                         # all validated scenarios
    uv run python -m kessler_sim.validation --scenario fengyun_1c_2007

Validated scenarios (see ``VALIDATIONS``):
* ``iridium_cosmos_2009``: Iridium 33 x Cosmos 2251 collision, real debris from the
  iridium-33-debris and cosmos-2251-debris TLE groups;
* ``fengyun_1c_2007``: anti-satellite test, real debris from the fengyun-1c-debris group
  (only Fengyun-1C's own debris: the interceptor's fragments fall back, in reality and
  in the simulation).

Real debris: the CelesTrak TLE groups in ``data/raw/tle/``, i.e. the fragments still in
orbit when the data was downloaded, 17-20 years after the events. Simulated debris:
right after the collision.

What can and cannot be compared:
* inclination barely changes over the years -> the strongest comparison;
* perigee / apogee and fragment counts are changed by 17 years of atmospheric drag,
  which the simulation does not model for fragments -> qualitative only.

The simulation is not adjusted to the results; differences are documented as
limitations. Outputs per scenario, with prefix e.g. ``validation_iridium_cosmos``:
    results/figures/<prefix>_inclination.png
    results/figures/<prefix>_gabbard.png
    results/figures/<prefix>_perigee_apogee.png
    results/<prefix>.json
    results/<prefix>.md
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from sgp4.api import Satrec

from . import orbits
from .catalogue import PROJECT_ROOT, jd_to_datetime
from .scenario import load_scenario

TLE_DIR = PROJECT_ROOT / "data" / "raw" / "tle"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

VALIDATIONS = {
    # scenario: title, output prefix, parents compared (NORAD ID, name, real-debris TLE file)
    "iridium_cosmos_2009": {
        "title": "Iridium 33 x Cosmos 2251",
        "prefix": "validation_iridium_cosmos",
        "parents": ((24946, "IRIDIUM 33", "iridium-33-debris.tle"),
                    (22675, "COSMOS 2251", "cosmos-2251-debris.tle")),
    },
    "fengyun_1c_2007": {
        "title": "Fengyun-1C anti-satellite test",
        "prefix": "validation_fengyun_1c",
        "parents": ((25730, "FENGYUN 1C", "fengyun-1c-debris.tle"),),
    },
}
PARENTS = VALIDATIONS["iridium_cosmos_2009"]["parents"]

# categorical slots 1 and 2 of the reference palette (validated: CVD and normal-vision pass on white)
SIM_COLOR = "#257de0"
REAL_COLOR = "#f47230"
TEXT = "#1a1a1a"
TEXT_SECONDARY = "#5c5c5c"
GRID = "#e6e6e6"


@dataclass
class Population:
    """Orbits of one set of fragments."""

    label: str
    inclination_deg: np.ndarray
    perigee_km: np.ndarray
    apogee_km: np.ndarray
    period_min: np.ndarray
    epoch: str = ""

    @property
    def count(self) -> int:
        return len(self.inclination_deg)


def _population(label, a_km, e, i_rad, epoch="") -> Population:
    a_km, e = np.asarray(a_km, float), np.asarray(e, float)
    return Population(
        label=label,
        inclination_deg=np.degrees(np.asarray(i_rad, float)),
        perigee_km=a_km * (1.0 - e) - orbits.R_MEAN,
        apogee_km=a_km * (1.0 + e) - orbits.R_MEAN,
        period_min=2.0 * np.pi * np.sqrt(a_km ** 3 / orbits.MU) / 60.0,
        epoch=epoch,
    )


def load_real_debris(tle_file: Path) -> Population:
    """Mean orbits of the real fragments from a TLE file (at each TLE's own epoch)."""
    lines = [ln.strip() for ln in Path(tle_file).read_text(encoding="utf-8").splitlines() if ln.strip()]
    sats = [Satrec.twoline2rv(lines[k + 1], lines[k + 2]) for k in range(len(lines) - 2)
            if lines[k + 1].startswith("1 ") and lines[k + 2].startswith("2 ")]
    if not sats:
        raise ValueError(f"no TLE records in {tle_file}")
    n = np.array([s.no_kozai for s in sats]) / 60.0                         # rad/s
    a = (orbits.MU / n ** 2) ** (1.0 / 3.0)
    epochs = sorted(jd_to_datetime(s.jdsatepoch + s.jdsatepochF) for s in sats)
    epoch = f"{epochs[0]:%Y-%m-%d} .. {epochs[-1]:%Y-%m-%d}"
    return _population("real", a, [s.ecco for s in sats], [s.inclo for s in sats], epoch)


def simulated_debris(sim, event, parent_id: int) -> Population:
    """Orbits of the simulated fragments of one parent, right after the collision."""
    ids = event.fragment_ids[event.fragment_origin_ids == parent_id]
    rows = [sim.table.row(i) for i in ids]
    el = sim.table.elements.take(rows)
    return _population("simulated", el.a, el.e, el.i, f"{sim.utc_at(event.t):%Y-%m-%d %H:%M} UTC")


def describe(p: Population) -> dict:
    def pct(x, q):
        return float(np.percentile(x, q))
    return {
        "count": p.count,
        "epoch": p.epoch,
        "inclination_mean_deg": float(p.inclination_deg.mean()),
        "inclination_std_deg": float(p.inclination_deg.std()),
        "perigee_p5_km": pct(p.perigee_km, 5), "perigee_median_km": pct(p.perigee_km, 50),
        "perigee_p95_km": pct(p.perigee_km, 95),
        "apogee_p5_km": pct(p.apogee_km, 5), "apogee_median_km": pct(p.apogee_km, 50),
        "apogee_p95_km": pct(p.apogee_km, 95),
        "apogee_above_1500km_fraction": float(np.mean(p.apogee_km > 1500.0)),
    }


def compare(tle_dir: Path = TLE_DIR, scenario: str = "iridium_cosmos_2009"):
    """Returns (summary dict, {parent name: (simulated, real)})."""
    spec = VALIDATIONS[scenario]
    sim = load_scenario(scenario)
    event = sim.events[0]
    populations, parents = {}, {}
    for pid, name, fname in spec["parents"]:
        simulated = simulated_debris(sim, event, pid)
        real = load_real_debris(Path(tle_dir) / fname)
        populations[name] = (simulated, real)
        s, r = describe(simulated), describe(real)
        parents[name] = {
            "norad_id": pid, "simulated": s, "real": r,
            "inclination_mean_difference_deg": s["inclination_mean_deg"] - r["inclination_mean_deg"],
            "inclination_std_ratio": s["inclination_std_deg"] / r["inclination_std_deg"],
            "apogee_p95_difference_km": s["apogee_p95_km"] - r["apogee_p95_km"],
        }
    compared = {pid for pid, _, _ in spec["parents"]}
    summary = {
        "scenario": scenario,
        "title": spec["title"],
        "collision_utc": f"{sim.utc_at(event.t):%Y-%m-%d %H:%M:%S}",
        "simulated_fragments": int(event.n_fragments),
        "simulated_reentered_immediately": int(event.reentered_count),
        "simulated_fragments_not_compared": {
            str(sim.table.names[sim.table.row(pid)]): int(np.sum(event.fragment_origin_ids == pid))
            for pid in (event.id_a, event.id_b) if pid not in compared
        },
        "parents": parents,
        "notes": [
            "Real debris is the population still in orbit at the TLE epochs (17-20 years after the event);"
            " simulated debris is the population right after the collision.",
            "Inclination is nearly unaffected by the years since the event: the strongest comparison.",
            "Perigee, apogee and counts are changed by atmospheric drag, which the simulation does not model"
            " for fragments: qualitative comparison only.",
        ],
    }
    return summary, populations


# -- plots --------------------------------------------------------------------------

def _style(ax):
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_SECONDARY)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9)
    ax.xaxis.label.set_color(TEXT)
    ax.yaxis.label.set_color(TEXT)
    ax.title.set_color(TEXT)


def _hist_pair(ax, sim_values, real_values, bins, xlabel, title):
    """Real = filled bars, simulated = outline (secondary encoding besides colour)."""
    ax.hist(real_values, bins=bins, density=True, color=REAL_COLOR, alpha=0.55, label="real debris",
            edgecolor="white", linewidth=0.6)
    ax.hist(sim_values, bins=bins, density=True, histtype="step", color=SIM_COLOR, linewidth=2.0,
            label="simulated")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("share of fragments (density)")
    ax.set_title(title, fontsize=11, loc="left")
    _style(ax)
    ax.legend(frameon=False, fontsize=9, labelcolor=TEXT)


def plot_inclination(populations, path: Path):
    import matplotlib.pyplot as plt
    n = len(populations)
    fig, axes = plt.subplots(1, n, figsize=(max(6 * n, 8), 4.4), squeeze=False)
    for ax, (name, (s, r)) in zip(axes[0], populations.items()):
        centre = r.inclination_deg.mean()
        both = np.concatenate([s.inclination_deg, r.inclination_deg]) - centre
        half = max(1.5, float(np.percentile(np.abs(both), 98)))       # covers 98 % of both, at least ±1.5°
        bins = np.linspace(centre - half, centre + half, 61)
        _hist_pair(ax, s.inclination_deg, r.inclination_deg, bins, "inclination [deg]",
                   f"{name}: inclination of fragments")
        ax.annotate(f"mean: simulated {s.inclination_deg.mean():.2f}°, real {r.inclination_deg.mean():.2f}°\n"
                    f"spread (std): simulated {s.inclination_deg.std():.2f}°, real {r.inclination_deg.std():.2f}°",
                    xy=(0.02, 0.97), xycoords="axes fraction", va="top", fontsize=9, color=TEXT_SECONDARY)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def plot_gabbard(populations, path: Path):
    import matplotlib.pyplot as plt
    n = len(populations)
    fig, axes = plt.subplots(1, n, figsize=(max(6 * n, 9), 5), sharey=True, squeeze=False)
    axes = axes[0]
    for ax, (name, (s, r)) in zip(axes, populations.items()):
        # simulated first, real on top (fewer points, filled)
        for pop, color, label, filled in ((s, SIM_COLOR, "simulated", False), (r, REAL_COLOR, "real debris", True)):
            face = color if filled else "none"
            ax.scatter(pop.period_min, pop.apogee_km, s=12, marker="^", facecolors=face, edgecolors=color,
                       linewidths=0.8, alpha=0.8, label=f"{label}: apogee")
            ax.scatter(pop.period_min, pop.perigee_km, s=12, marker="v", facecolors=face, edgecolors=color,
                       linewidths=0.8, alpha=0.8, label=f"{label}: perigee")
        ax.set_xlim(92, 125)
        ax.set_ylim(0, 3000)
        ax.set_xlabel("orbital period [min]")
        ax.set_title(f"{name}: Gabbard diagram", fontsize=11, loc="left")
        _style(ax)
        ax.legend(frameon=False, fontsize=8, labelcolor=TEXT, loc="upper left")
    axes[0].set_ylabel("altitude [km]")
    fig.text(0.5, -0.04, "Simulated: right after the collision. Real: still in orbit 17-20 years later;\n"
             "drag has removed low orbits and lowered apogees. Axes cut at 125 min / 3000 km.",
             ha="center", fontsize=9, color=TEXT_SECONDARY)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_perigee_apogee(populations, path: Path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(populations), 2, figsize=(12, 4 * len(populations)), squeeze=False)
    perigee_bins = np.linspace(0, 1200, 49)
    apogee_bins = np.linspace(0, 3000, 61)
    for row, (name, (s, r)) in zip(axes, populations.items()):
        _hist_pair(row[0], np.clip(s.perigee_km, 0, 1200), np.clip(r.perigee_km, 0, 1200), perigee_bins,
                   "perigee altitude [km]", f"{name}: perigee")
        _hist_pair(row[1], np.clip(s.apogee_km, 0, 3000), np.clip(r.apogee_km, 0, 3000), apogee_bins,
                   "apogee altitude [km]", f"{name}: apogee")
    fig.text(0.5, -0.01, "Apogees above 3000 km are shown in the last bin. Real debris: 17-20 years after the "
             "collision (drag has removed low orbits); simulated: right after it.", ha="center", fontsize=9,
             color=TEXT_SECONDARY)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


# -- report -------------------------------------------------------------------------

def markdown_table(summary: dict) -> str:
    lines = [
        f"# Validation: {summary['title']}",
        "",
        f"Simulated collision: {summary['collision_utc']} UTC, {summary['simulated_fragments']} fragments "
        f"(≥ 10 cm) in orbit, {summary['simulated_reentered_immediately']} re-entered immediately.",
    ]
    if summary.get("simulated_fragments_not_compared"):
        lines.append("Not compared (no real data): " + ", ".join(
            f"{n}: {k} fragments in orbit" for n, k in summary["simulated_fragments_not_compared"].items()) + ".")
    lines += [
        "",
        "| Quantity | " + " | ".join(f"{n}: simulated | {n}: real" for n in summary["parents"]) + " |",
        "|---|" + "---|---|" * len(summary["parents"]),
    ]
    rows = [
        ("Fragments", "count", "{:.0f}"),
        ("Inclination mean [deg]", "inclination_mean_deg", "{:.2f}"),
        ("Inclination std [deg]", "inclination_std_deg", "{:.2f}"),
        ("Perigee 5–95 % [km]", ("perigee_p5_km", "perigee_p95_km"), "{:.0f}–{:.0f}"),
        ("Apogee 5–95 % [km]", ("apogee_p5_km", "apogee_p95_km"), "{:.0f}–{:.0f}"),
        ("Apogee above 1500 km", "apogee_above_1500km_fraction", "{:.0%}"),
    ]
    for title, key, fmt in rows:
        cells = []
        for p in summary["parents"].values():
            for side in ("simulated", "real"):
                d = p[side]
                cells.append(fmt.format(*(d[k] for k in key)) if isinstance(key, tuple) else fmt.format(d[key]))
        lines.append(f"| {title} | " + " | ".join(cells) + " |")
    lines += ["", "Real data epochs: " + ", ".join(f"{n}: {p['real']['epoch']}" for n, p in summary["parents"].items()),
              "", "Notes:"] + [f"- {n}" for n in summary["notes"]]
    return "\n".join(lines) + "\n"


def run(tle_dir: Path = TLE_DIR, figures_dir: Path = FIGURES_DIR, results_dir: Path = RESULTS_DIR,
        scenario: str = "iridium_cosmos_2009") -> dict:
    import matplotlib
    matplotlib.use("Agg")
    summary, populations = compare(tle_dir, scenario)
    figures_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    prefix = VALIDATIONS[scenario]["prefix"]
    outputs = {
        "inclination": figures_dir / f"{prefix}_inclination.png",
        "gabbard": figures_dir / f"{prefix}_gabbard.png",
        "perigee_apogee": figures_dir / f"{prefix}_perigee_apogee.png",
        "json": results_dir / f"{prefix}.json",
        "markdown": results_dir / f"{prefix}.md",
    }
    plot_inclination(populations, outputs["inclination"])
    plot_gabbard(populations, outputs["gabbard"])
    plot_perigee_apogee(populations, outputs["perigee_apogee"])
    outputs["json"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    outputs["markdown"].write_text(markdown_table(summary), encoding="utf-8")
    summary["outputs"] = {k: str(v) for k, v in outputs.items()}
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description="Compare simulated debris with the real debris of the same event")
    ap.add_argument("--scenario", choices=[*VALIDATIONS, "all"], default="all")
    ap.add_argument("--figures-dir", type=Path, default=FIGURES_DIR)
    ap.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    args = ap.parse_args(argv)
    for scenario in (VALIDATIONS if args.scenario == "all" else [args.scenario]):
        summary = run(figures_dir=args.figures_dir, results_dir=args.results_dir, scenario=scenario)
        print(markdown_table(summary))
        for path in summary["outputs"].values():
            print("written:", path)
        print()


if __name__ == "__main__":
    main()
