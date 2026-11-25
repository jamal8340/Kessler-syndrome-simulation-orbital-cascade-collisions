import json

import numpy as np
import pytest

from kessler_sim import validation as val


@pytest.fixture(scope="module")
def comparison():
    return val.compare()


def _records_in(path):
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    return sum(1 for k in range(len(lines) - 2) if lines[k + 1].startswith("1 ") and lines[k + 2].startswith("2 "))


ALL_PARENTS = [p for spec in val.VALIDATIONS.values() for p in spec["parents"]]


@pytest.mark.parametrize("pid,name,fname", ALL_PARENTS)
def test_real_debris_matches_tle_file(pid, name, fname):
    path = val.TLE_DIR / fname
    pop = val.load_real_debris(path)
    assert pop.count == _records_in(path)
    assert np.all(pop.apogee_km >= pop.perigee_km)
    assert np.all((pop.period_min > 85) & (pop.period_min < 200))


def test_population_from_known_orbit():
    # circular orbit at 789 km, 86.4 deg: perigee = apogee = 789 km, period ~100.4 min
    pop = val._population("x", [val.orbits.R_MEAN + 789.0], [0.0], [np.radians(86.4)])
    assert pop.perigee_km[0] == pytest.approx(789.0)
    assert pop.apogee_km[0] == pytest.approx(789.0)
    assert pop.inclination_deg[0] == pytest.approx(86.4)
    assert pop.period_min[0] == pytest.approx(100.4, abs=0.2)


def test_simulated_fragments_split_by_parent(comparison):
    summary, populations = comparison
    counts = [p["simulated"]["count"] for p in summary["parents"].values()]
    assert sum(counts) == summary["simulated_fragments"]
    assert set(populations) == {"IRIDIUM 33", "COSMOS 2251"}


def test_inclination_means_agree_with_real_debris(comparison):
    """The strongest check: fragments keep their parent's orbital plane."""
    summary, _ = comparison
    for p in summary["parents"].values():
        assert abs(p["inclination_mean_difference_deg"]) < 0.2


def test_documented_limitation_wider_spread(comparison):
    """Known, accepted limitation: simulated kicks are too strong, so the spread is wider than real."""
    summary, _ = comparison
    for p in summary["parents"].values():
        assert p["inclination_std_ratio"] > 1.0
        assert p["simulated"]["apogee_p95_km"] > p["real"]["apogee_p95_km"]


@pytest.mark.parametrize("scenario", list(val.VALIDATIONS))
def test_run_writes_all_outputs(tmp_path, scenario):
    summary = val.run(figures_dir=tmp_path, results_dir=tmp_path, scenario=scenario)
    prefix = val.VALIDATIONS[scenario]["prefix"]
    for key in ("inclination", "gabbard", "perigee_apogee", "json", "markdown"):
        path = tmp_path / summary["outputs"][key].split("/")[-1]
        assert path.name.startswith(prefix)
        assert path.is_file() and path.stat().st_size > 0
    data = json.loads((tmp_path / f"{prefix}.json").read_text())
    assert set(data["parents"]) == {name for _, name, _ in val.VALIDATIONS[scenario]["parents"]}
    assert "| Inclination mean [deg] |" in (tmp_path / f"{prefix}.md").read_text()


def test_fengyun_comparison():
    summary, _ = val.compare(scenario="fengyun_1c_2007")
    fy = summary["parents"]["FENGYUN 1C"]
    assert abs(fy["inclination_mean_difference_deg"]) < 0.5        # looser than Iridium x Cosmos (see VALIDATION.md)
    assert fy["inclination_std_ratio"] > 1.0                       # same accepted limitation: kicks too strong
    assert summary["simulated_fragments_not_compared"] == {"ASAT INTERCEPTOR": 0}   # all fell back
