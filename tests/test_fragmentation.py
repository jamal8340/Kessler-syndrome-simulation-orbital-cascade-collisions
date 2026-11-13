import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "collision_physics"))
import fragmentation as fr

# Iridium 33 (~560 kg) i Cosmos 2251 (~900 kg), v_rel ~ 11.7 km/s, 10.02.2009
COSMOS = fr.Body(900.0, 1.5, (7.0, 0.0, 0.0))
IRIDIUM = fr.Body(560.0, 1.2, (-3.0, 8.9, 0.0))


def run(seed=1, **kw):
    return fr.simulate_collision(COSMOS, IRIDIUM, seed=seed, **kw)


def test_emr_known_values():
    e, emr, cat = fr.collision_energetics(100.0, 1.0, 10.0)   # 1 kg w 100 kg
    assert abs(e - 5e7) < 1e-6 and abs(emr - 500.0) < 1e-9 and cat
    e, emr, cat = fr.collision_energetics(1000.0, 0.01, 10.0)  # 10 g w 1000 kg
    assert abs(emr - 0.5) < 1e-9 and not cat


def test_catastrophic_and_count_formula():
    r = run()
    assert r.catastrophic
    expected = 0.1 * 1460.0 ** 0.75 * 0.1 ** -1.71
    assert abs(r.n_expected - round(expected)) <= 1


def test_mass_budget_not_exceeded():
    r = run()
    assert r.mass_kg.sum() <= r.mass_involved_kg + 1e-9
    assert abs(r.mass_kg.sum() + r.remnant_mass_kg - 1460.0) < 1e-6


def test_momentum_conserved():
    r = run()
    p_before = 900.0 * np.array(COSMOS.velocity_kms) + 560.0 * np.array(IRIDIUM.velocity_kms)
    p_after = (r.mass_kg[:, None] * r.velocity_kms).sum(axis=0) + r.remnant_mass_kg * r.remnant_velocity_kms
    assert np.allclose(p_before, p_after, rtol=1e-10, atol=1e-8)


def test_all_fragments_above_lmin():
    r = run()
    assert r.lc_m.min() >= 0.10 - 1e-12


def test_size_distribution_slope():
    r = run(l_min_m=0.01)               # więcej fragmentów dla stabilnego dopasowania
    grid = np.logspace(-1.8, -0.7, 15)
    counts = np.array([(r.lc_m >= g).sum() for g in grid])
    slope = np.polyfit(np.log10(grid), np.log10(counts), 1)[0]
    assert abs(slope + 1.71) < 0.15, slope


def test_reproducible_with_seed():
    a, b = run(seed=7), run(seed=7)
    assert np.array_equal(a.mass_kg, b.mass_kg)
    assert not np.array_equal(a.mass_kg, run(seed=8).mass_kg)


def test_energy_not_exceeded():
    assert run().dv_kinetic_fraction < 1.0


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("[PASS]", name)
    r = run()
    print(f"\nIridium/Cosmos: EMR={r.emr_jg:.0f} J/g, N(wzór)={r.n_expected}, "
          f"zachowane={r.n_fragments}, masa fragmentów={r.mass_kg.sum():.0f} kg, "
          f"nieśledzona={r.remnant_mass_kg:.0f} kg, Δv/E_CM={r.dv_kinetic_fraction:.1e}")