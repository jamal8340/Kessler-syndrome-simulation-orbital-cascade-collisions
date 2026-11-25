import sys
from pathlib import Path

import numpy as np
import pytest

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
    assert r.mass_kg.sum() <= r.mass_budget_kg + 1e-9
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


# --- przypadki niekatastroficzne i brzegowe --------------------------------

TARGET = fr.Body(1000.0, 1.5, (7.0, 0.0, 0.0))
SMALL_PROJECTILE = fr.Body(0.01, 0.05, (0.0, 7.0, 0.0))   # 10 g, v_rel ~ 9.9 km/s


def _momentum(r):
    return (r.mass_kg[:, None] * r.velocity_kms).sum(axis=0) + r.remnant_mass_kg * r.remnant_velocity_kms


def test_non_catastrophic_mass_involved_and_budget():
    r = fr.simulate_collision(TARGET, SMALL_PROJECTILE, seed=1)
    v_rel = np.linalg.norm(np.subtract(TARGET.velocity_kms, SMALL_PROJECTILE.velocity_kms))
    assert not r.catastrophic
    assert r.mass_involved_kg == pytest.approx(0.01 * v_rel ** 2)        # M = m_p * v^2 [km/s]
    assert r.mass_budget_kg == pytest.approx(min(r.mass_involved_kg, 1000.0))
    assert r.mass_kg.sum() <= r.mass_budget_kg + 1e-9
    assert r.n_expected == round(0.1 * r.mass_involved_kg ** 0.75 * 0.1 ** -1.71)


def test_non_catastrophic_momentum_conserved():
    r = fr.simulate_collision(TARGET, SMALL_PROJECTILE, seed=2)
    p_before = 1000.0 * np.array(TARGET.velocity_kms) + 0.01 * np.array(SMALL_PROJECTILE.velocity_kms)
    assert np.allclose(_momentum(r), p_before, rtol=1e-10, atol=1e-10)


def test_total_mass_conserved_non_catastrophic():
    r = fr.simulate_collision(TARGET, SMALL_PROJECTILE, seed=3)
    assert r.mass_kg.sum() + r.remnant_mass_kg == pytest.approx(1000.01)


def test_result_independent_of_argument_order():
    a, b = run(seed=5), fr.simulate_collision(IRIDIUM, COSMOS, seed=5)
    assert np.array_equal(a.mass_kg, b.mass_kg)
    assert np.array_equal(a.velocity_kms, b.velocity_kms)


def test_bodies_too_small_for_lmin_give_no_fragments():
    # l_max = 2 * max(radius) = 0.08 m < L_MIN = 0.1 m
    a = fr.Body(1.0, 0.04, (7.0, 0.0, 0.0))
    b = fr.Body(1.0, 0.04, (0.0, 7.0, 0.0))
    r = fr.simulate_collision(a, b, seed=1)
    assert r.n_fragments == 0
    assert r.remnant_mass_kg == pytest.approx(2.0)
    assert np.allclose(r.remnant_velocity_kms, [3.5, 3.5, 0.0])
    assert r.dv_kinetic_fraction == 0.0


def test_zero_relative_velocity_gives_no_fragments():
    a = fr.Body(500.0, 1.5, (7.0, 0.0, 0.0))
    b = fr.Body(100.0, 1.0, (7.0, 0.0, 0.0))
    r = fr.simulate_collision(a, b, seed=1)
    assert not r.catastrophic and r.n_fragments == 0
    assert r.remnant_mass_kg == pytest.approx(600.0)


def test_fragment_properties_physical():
    r = run()
    assert np.all(r.mass_kg > 0)
    assert np.all(r.am_m2_per_kg > 0)
    assert np.all(np.isfinite(r.velocity_kms))
    assert r.lc_m.max() <= 2.0 * max(COSMOS.radius_m, IRIDIUM.radius_m) + 1e-12
    assert r.velocity_kms.shape == (r.n_fragments, 3)


def test_fragment_area_matches_size_and_am():
    r = run()
    # A = 0.540424 * Lc^2 oraz m = A / (A/M)
    assert np.allclose(r.mass_kg * r.am_m2_per_kg, 0.540424 * r.lc_m ** 2)


def test_fragments_move_around_centre_of_mass():
    r = run()
    v_cm = (900.0 * np.array(COSMOS.velocity_kms) + 560.0 * np.array(IRIDIUM.velocity_kms)) / 1460.0
    assert np.allclose(r.v_cm_kms, v_cm)
    dv = r.velocity_kms - r.v_cm_kms
    assert np.allclose((r.mass_kg[:, None] * dv).sum(axis=0), 0.0, atol=1e-9)


def test_delta_v_distribution_follows_am_relation():
    # log10(dv [m/s]) ~ N(0.9*chi + 2.9, 0.4): reszty po odjęciu średniej mają sigma ~ 0.4
    a = fr.Body(1000.0, 2.0, (7.0, 0.0, 0.0))
    b = fr.Body(1000.0, 2.0, (-7.0, 0.0, 0.0))
    rng = np.random.default_rng(0)
    chi = rng.uniform(-2, 0, 20000)
    dv = np.linalg.norm(fr.sample_delta_v_kms(chi, rng), axis=1) * 1000.0
    resid = np.log10(dv) - (0.9 * chi + 2.9)
    assert abs(resid.mean()) < 0.02
    assert abs(resid.std() - 0.4) < 0.02
    del a, b


def test_delta_v_directions_isotropic():
    rng = np.random.default_rng(1)
    d = fr.sample_delta_v_kms(np.zeros(20000), rng)
    u = d / np.linalg.norm(d, axis=1, keepdims=True)
    assert np.all(np.abs(u.mean(axis=0)) < 0.03)


def test_sample_lengths_within_bounds():
    rng = np.random.default_rng(2)
    lc = fr.sample_lengths(10000, 0.1, 3.0, rng)
    assert lc.min() >= 0.1 and lc.max() <= 3.0


def test_bigger_impact_gives_more_fragments():
    slow = fr.simulate_collision(TARGET, fr.Body(0.5, 0.1, (6.0, 0.0, 0.0)), seed=1)    # 1 km/s
    fast = fr.simulate_collision(TARGET, fr.Body(0.5, 0.1, (-3.0, 0.0, 0.0)), seed=1)   # 10 km/s
    assert fast.n_expected > slow.n_expected


@pytest.mark.xfail(strict=True, reason=(
    "Znany problem modelu: rozkład dv w NASA SBM nie zależy od energii zderzenia, więc dla "
    "słabych zderzeń niekatastroficznych energia kinetyczna fragmentów może przekroczyć "
    "energię dostępną w układzie środka masy."))
def test_energy_not_exceeded_for_weak_non_catastrophic_impact():
    # 4 g w 15 kg przy 7 km/s; występuje dla ~3 na 2000 ziaren (tu seed=18 -> ~1.6x)
    target = fr.Body(15.0, 0.5, (7.0, 0.0, 0.0))
    projectile = fr.Body(0.004, 0.06, (0.0, 0.0, 0.0))
    r = fr.simulate_collision(target, projectile, seed=18)
    assert r.dv_kinetic_fraction < 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
    r = run()
    print(f"\nIridium/Cosmos: EMR={r.emr_jg:.0f} J/g, N(wzór)={r.n_expected}, "
          f"zachowane={r.n_fragments}, masa fragmentów={r.mass_kg.sum():.0f} kg, "
          f"nieśledzona={r.remnant_mass_kg:.0f} kg, Δv/E_CM={r.dv_kinetic_fraction:.1e}")