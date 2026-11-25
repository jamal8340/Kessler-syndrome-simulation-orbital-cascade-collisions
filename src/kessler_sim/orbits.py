"""Analytic orbit model: two-body motion plus secular J2 drift.

Every object is stored as a set of mean orbital elements referenced to a time ``t_ref``.
The position at any time t is computed directly (no time-stepping):

    RAAN(t)  = RAAN0  + dRAAN/dt * (t - t_ref)
    argp(t)  = argp0  + dargp/dt * (t - t_ref)
    M(t)     = M0     + dM/dt    * (t - t_ref) + 1/2 * ndot * (t - t_ref)^2

with the classic first-order J2 secular rates. ``ndot`` is an optional per-object
mean-motion decay rate (taken from SGP4 for catalogue objects, zero for fragments): it
is the only drag effect in the model. Short-period J2 terms, higher harmonics and
lunar/solar perturbations are ignored (documented simplification).

Units: km, s, rad. Constants are WGS-72, the same as SGP4, so mean elements taken
from SGP4 are used consistently.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MU = 398600.8               # km^3/s^2 (WGS-72, as in SGP4)
R_EQ = 6378.135             # km, equatorial radius (WGS-72)
J2 = 0.001082616
R_MEAN = 6371.0             # km, mean radius used for reporting altitudes

TWO_PI = 2.0 * np.pi
_EPS = 1e-10


@dataclass
class Elements:
    """Arrays of mean elements (all shape (n,))."""

    a: np.ndarray       # semi-major axis [km]
    e: np.ndarray       # eccentricity
    i: np.ndarray       # inclination [rad]
    raan: np.ndarray    # right ascension of ascending node [rad]
    argp: np.ndarray    # argument of perigee [rad]
    m: np.ndarray       # mean anomaly at t_ref [rad]
    ndot: np.ndarray | None = None   # mean-motion decay [rad/s^2]; None = no decay

    FIELDS = ("a", "e", "i", "raan", "argp", "m", "ndot")

    def __post_init__(self):
        if self.ndot is None:
            self.ndot = np.zeros_like(np.asarray(self.a, dtype=float))

    def __len__(self):
        return len(self.a)

    @staticmethod
    def empty() -> "Elements":
        return Elements(*(np.empty(0) for _ in Elements.FIELDS))

    @staticmethod
    def concat(parts: list["Elements"]) -> "Elements":
        return Elements(*(np.concatenate([getattr(p, f) for p in parts]) for f in Elements.FIELDS))

    def take(self, idx) -> "Elements":
        return Elements(*(getattr(self, f)[idx] for f in Elements.FIELDS))


def _as2d(x):
    x = np.asarray(x, dtype=float)
    return x[None, :] if x.ndim == 1 else x


def secular_rates(el: Elements):
    """First-order J2 secular rates (dRAAN/dt, dargp/dt, dM/dt) in rad/s."""
    n = np.sqrt(MU / el.a ** 3)
    p = el.a * (1.0 - el.e ** 2)
    k = 1.5 * J2 * (R_EQ / p) ** 2 * n
    cos_i = np.cos(el.i)
    raan_dot = -k * cos_i
    argp_dot = 0.5 * k * (5.0 * cos_i ** 2 - 1.0)
    m_dot = n + 0.5 * k * np.sqrt(1.0 - el.e ** 2) * (3.0 * cos_i ** 2 - 1.0)
    return raan_dot, argp_dot, m_dot


def solve_kepler(m, e, iterations=12):
    """Eccentric anomaly E from mean anomaly M (vectorised Newton iteration)."""
    m = np.mod(m, TWO_PI)
    ecc_anom = np.where(e < 0.8, m, np.pi)
    for _ in range(iterations):
        f = ecc_anom - e * np.sin(ecc_anom) - m
        ecc_anom = ecc_anom - f / (1.0 - e * np.cos(ecc_anom))
    return ecc_anom


def _perifocal_to_inertial(raan, i, argp):
    """Columns P and Q of the perifocal->inertial rotation matrix, shape (n, 3) each."""
    cO, sO = np.cos(raan), np.sin(raan)
    ci, si = np.cos(i), np.sin(i)
    cw, sw = np.cos(argp), np.sin(argp)
    p_vec = np.stack([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si], axis=-1)
    q_vec = np.stack([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si], axis=-1)
    return p_vec, q_vec


def elements_to_state(el: Elements, dt, with_j2=True):
    """Position and velocity (n, 3) after time ``dt`` [s] (scalar or (n,)) from t_ref."""
    dt = np.broadcast_to(np.asarray(dt, dtype=float), el.a.shape)
    if with_j2:
        raan_dot, argp_dot, m_dot = secular_rates(el)
    else:
        raan_dot = argp_dot = 0.0
        m_dot = np.sqrt(MU / el.a ** 3)
    raan = el.raan + raan_dot * dt
    argp = el.argp + argp_dot * dt
    m = el.m + m_dot * dt + 0.5 * el.ndot * dt ** 2

    ecc_anom = solve_kepler(m, el.e)
    cos_e, sin_e = np.cos(ecc_anom), np.sin(ecc_anom)
    root = np.sqrt(1.0 - el.e ** 2)
    x_pf = el.a * (cos_e - el.e)
    y_pf = el.a * root * sin_e
    n = np.sqrt(MU / el.a ** 3)
    denom = 1.0 - el.e * cos_e
    vx_pf = -el.a * n * sin_e / denom
    vy_pf = el.a * n * root * cos_e / denom

    p_vec, q_vec = _perifocal_to_inertial(raan, el.i, argp)
    r = x_pf[:, None] * p_vec + y_pf[:, None] * q_vec
    v = vx_pf[:, None] * p_vec + vy_pf[:, None] * q_vec
    return r, v


def state_to_elements(r, v):
    """Osculating elements from position/velocity arrays (n, 3).

    Singular cases are handled with fixed conventions: for (near-)circular orbits
    argp = 0 and the mean anomaly is the argument of latitude; for (near-)equatorial
    orbits RAAN = 0 and angles are measured from the x axis.

    Returns (Elements, bound) where ``bound`` is False for parabolic/hyperbolic states.
    """
    r = _as2d(r)
    v = _as2d(v)
    r_mag = np.linalg.norm(r, axis=1)
    v_sq = np.einsum("ij,ij->i", v, v)
    h = np.cross(r, v)
    h_mag = np.linalg.norm(h, axis=1)
    h_hat = h / h_mag[:, None]

    energy = 0.5 * v_sq - MU / r_mag
    bound = energy < 0
    a = np.where(bound, -MU / (2.0 * np.where(bound, energy, -1.0)), np.inf)

    e_vec = ((v_sq - MU / r_mag)[:, None] * r - np.einsum("ij,ij->i", r, v)[:, None] * v) / MU
    e = np.linalg.norm(e_vec, axis=1)

    i = np.arccos(np.clip(h_hat[:, 2], -1.0, 1.0))

    node = np.stack([-h[:, 1], h[:, 0], np.zeros(len(h))], axis=1)   # z x h
    node_mag = np.linalg.norm(node, axis=1)
    equatorial = node_mag < _EPS * h_mag
    node_hat = np.where(equatorial[:, None], np.array([1.0, 0.0, 0.0]),
                        node / np.where(equatorial, 1.0, node_mag)[:, None])
    raan = np.where(equatorial, 0.0, np.arctan2(node_hat[:, 1], node_hat[:, 0]))

    perp = np.cross(h_hat, node_hat)                 # in-plane, 90 deg ahead of the node
    u = np.arctan2(np.einsum("ij,ij->i", r, perp), np.einsum("ij,ij->i", r, node_hat))

    circular = e < 1e-9
    argp = np.where(circular, 0.0,
                    np.arctan2(np.einsum("ij,ij->i", e_vec, perp), np.einsum("ij,ij->i", e_vec, node_hat)))
    nu = u - argp

    e_safe = np.minimum(e, 0.999999999)
    ecc_anom = np.arctan2(np.sqrt(1.0 - e_safe ** 2) * np.sin(nu), e_safe + np.cos(nu))
    m = ecc_anom - e_safe * np.sin(ecc_anom)

    el = Elements(a, e, i, np.mod(raan, TWO_PI), np.mod(argp, TWO_PI), np.mod(m, TWO_PI))
    return el, bound & (e < 1.0)


def perigee_radius(el: Elements):
    return el.a * (1.0 - el.e)


def apogee_radius(el: Elements):
    return el.a * (1.0 + el.e)


def period(el: Elements):
    return TWO_PI * np.sqrt(el.a ** 3 / MU)
