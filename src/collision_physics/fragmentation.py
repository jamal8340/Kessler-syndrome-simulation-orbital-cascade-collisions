from __future__ import annotations
from dataclasses import dataclass
import numpy as np

# constants
CATASTROPHIC_EMR_JG = 40.0          # próg katastroficzności [J/g]
L_MIN_DEFAULT_M = 0.10              # dolna granica generowanych fragmentów [m]
LC_AM_VALIDITY_M = 0.11             # dolna granica ważności bimodalnego A/M [m]

N_COEFF = 0.1                       # N(>Lc) = 0.1 * M^0.75 * Lc^-1.71
N_MASS_EXP = 0.75
N_SIZE_EXP = 1.71

AREA_COEFF = 0.540424               # A = 0.540424 * Lc^2   (Lc >= 1.67 mm)

# delta v: log10(delta v[m/s]) ~ N(0.9*chi + 2.9, 0.4),  chi = log10(A/M)   (kolizja)
DV_MU_SLOPE, DV_MU_OFFSET, DV_SIGMA = 0.9, 2.9, 0.4

# A/M dla statku kosmicznego: kawałkami liniowe w lambda = log10(Lc)
# (wartości na krańcach przedziałów, np.interp trzyma stałe wartości poza nimi)
_AM_ALPHA = ([-1.95, 0.55], [0.0, 1.0])
_AM_MU1 = ([-1.1, 0.0], [-0.6, -0.95])
_AM_SIG1 = ([-1.3, -0.3], [0.1, 0.3])
_AM_MU2 = ([-0.7, -0.1], [-1.2, -2.0])
_AM_SIG2 = ([-0.5, -0.3], [0.5, 0.3])

@dataclass
class Body:
    mass_kg: float
    radius_m: float
    velocity_kms: tuple[float, float, float]


@dataclass
class BreakupResult:
    catastrophic: bool
    energy_j: float
    emr_jg: float
    mass_involved_kg: float          # M we wzorze na liczbę fragmentów
    mass_budget_kg: float            # limit sumy mas fragmentów
    n_expected: int                  # N(>= L_MIN) ze wzoru, przed przycięciem budżetem
    lc_m: np.ndarray                 # (n,)
    am_m2_per_kg: np.ndarray         # (n,)
    mass_kg: np.ndarray              # (n,)
    velocity_kms: np.ndarray         # (n, 3)
    remnant_mass_kg: float           # masa nieśledzona (pozostałość + fragmenty < L_MIN)
    remnant_velocity_kms: np.ndarray
    v_cm_kms: np.ndarray
    dv_kinetic_fraction: float       # energia kinetyczna Δv / energia dostępna w układzie CM

    @property
    def n_fragments(self) -> int:
        return len(self.mass_kg)


def collision_energetics(m_a_kg, m_b_kg, v_rel_kms):
    m_p, m_t = sorted((float(m_a_kg), float(m_b_kg)))
    if m_p <= 0:
        return 0.0, 0.0, False
    energy_j = 0.5 * m_p * (v_rel_kms * 1000.0) ** 2
    emr = energy_j / (m_t * 1000.0)
    return energy_j, emr, emr >= CATASTROPHIC_EMR_JG


def mass_involved(m_t_kg, m_p_kg, v_rel_kms, catastrophic):
    if catastrophic:
        return m_t_kg + m_p_kg
    return m_p_kg * v_rel_kms ** 2


def expected_fragment_count(m_involved_kg, l_min_m):
    return N_COEFF * m_involved_kg ** N_MASS_EXP * l_min_m ** (-N_SIZE_EXP)


def sample_lengths(n, l_min, l_max, rng):
    b = N_SIZE_EXP
    lo, hi = l_max ** (-b), l_min ** (-b)
    return (lo + rng.random(n) * (hi - lo)) ** (-1.0 / b)


def sample_log_area_to_mass(lc, rng):
    lam = np.log10(np.maximum(lc, LC_AM_VALIDITY_M))
    alpha = np.interp(lam, *_AM_ALPHA)
    mu1, s1 = np.interp(lam, *_AM_MU1), np.interp(lam, *_AM_SIG1)
    mu2, s2 = np.interp(lam, *_AM_MU2), np.interp(lam, *_AM_SIG2)
    first_mode = rng.random(len(lc)) < alpha
    return np.where(first_mode, rng.normal(mu1, s1), rng.normal(mu2, s2))


def sample_delta_v_kms(chi, rng):
    mu = DV_MU_SLOPE * chi + DV_MU_OFFSET
    speed_kms = 10.0 ** rng.normal(mu, DV_SIGMA) / 1000.0
    d = rng.normal(size=(len(chi), 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return d * speed_kms[:, None]


def simulate_collision(a: Body, b: Body, l_min_m=L_MIN_DEFAULT_M, seed=None) -> BreakupResult:
    rng = np.random.default_rng(seed)

    target, projectile = (a, b) if a.mass_kg >= b.mass_kg else (b, a)
    v_t = np.asarray(target.velocity_kms, float)
    v_p = np.asarray(projectile.velocity_kms, float)
    v_rel = float(np.linalg.norm(v_t - v_p))
    m_tot = target.mass_kg + projectile.mass_kg
    v_cm = (target.mass_kg * v_t + projectile.mass_kg * v_p) / m_tot

    energy_j, emr, cata = collision_energetics(target.mass_kg, projectile.mass_kg, v_rel)
    m_inv = mass_involved(target.mass_kg, projectile.mass_kg, v_rel, cata)
    # budżet masy: katastrofa -> cała masa; inaczej nie więcej niż masa celu
    budget = m_tot if cata else min(m_inv, target.mass_kg)

    n_exp = expected_fragment_count(m_inv, l_min_m)
    n = int(round(n_exp))
    l_max = 2.0 * max(target.radius_m, projectile.radius_m)

    empty = np.empty(0)
    if n == 0 or l_max <= l_min_m:
        return BreakupResult(cata, energy_j, emr, m_inv, budget, n, empty, empty, empty,
                             np.empty((0, 3)), m_tot, v_cm.copy(), v_cm, 0.0)

    lc = sample_lengths(n, l_min_m, l_max, rng)
    chi = sample_log_area_to_mass(lc, rng)
    area = AREA_COEFF * lc ** 2
    mass = area / 10.0 ** chi

    # budżet masy: gdy suma przekracza masę, usuwane są NAJWIĘKSZE fragmenty
    # (zostaje najdłuższy prefiks od najmniejszych, którego suma mieści się w budżecie)
    order = np.argsort(mass)
    keep = order[np.cumsum(mass[order]) <= budget]
    lc, chi, mass = lc[keep], chi[keep], mass[keep]

    dv = sample_delta_v_kms(chi, rng)
    if len(mass):
        dv -= (mass[:, None] * dv).sum(axis=0) / mass.sum()   # zachowanie pędu
    velocity = v_cm + dv

    ke_dv = 0.5 * float(np.sum(mass * (np.linalg.norm(dv, axis=1) * 1000.0) ** 2))
    mu_red = target.mass_kg * projectile.mass_kg / m_tot
    ke_cm = 0.5 * mu_red * (v_rel * 1000.0) ** 2

    return BreakupResult(
        catastrophic=cata, energy_j=energy_j, emr_jg=emr,
        mass_involved_kg=m_inv, mass_budget_kg=budget, n_expected=n,
        lc_m=lc, am_m2_per_kg=10.0 ** chi, mass_kg=mass, velocity_kms=velocity,
        remnant_mass_kg=m_tot - float(mass.sum()), remnant_velocity_kms=v_cm.copy(),
        v_cm_kms=v_cm, dv_kinetic_fraction=ke_dv / ke_cm if ke_cm > 0 else 0.0,
    )