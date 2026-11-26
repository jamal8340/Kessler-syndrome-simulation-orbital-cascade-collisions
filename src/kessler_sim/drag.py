"""Deterministic, orbit-averaged atmospheric drag for long-term LEO experiments.

Static, spherical, non-rotating atmosphere; no solar activity or explicit manoeuvres.
Optional ideal station keeping cancels drag on intact active satellites.
Units: km, s, rad for orbital elements, kg/m^3 for density, m^2/kg for A/M.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from . import orbits

DAY = 86400.0
# Nominal exponential atmosphere, Vallado-style reference density table.
# Log-linear interpolation makes density continuous between tabulated heights.
_HEIGHT = np.array([150, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000.])
_DENSITY = np.array([2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11,
                     9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13,
                     3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15])
_LOG_DENSITY = np.log(_DENSITY)
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(48)
_COS_E = np.cos((_NODES + 1) * np.pi / 2)
_WEIGHTS = _WEIGHTS / 2  # mean over E in [0, pi], symmetry over the other half


@dataclass(frozen=True)
class DragConfig:
    cd: float = 2.2
    density_scale: float = 1.0
    max_step_s: float = 90 * DAY
    max_axis_change_km: float = 5.0
    max_perigee_altitude_km: float = 2000.0
    reentry_altitude_km: float = 150.0
    stationkeeping: bool = True

    def __post_init__(self):
        if not all(np.isfinite(x) and x > 0 for x in
                   (self.cd, self.density_scale, self.max_step_s, self.max_axis_change_km,
                    self.max_perigee_altitude_km, self.reentry_altitude_km)):
            raise ValueError("drag parameters must be finite and positive")
        if not 150 <= self.reentry_altitude_km < self.max_perigee_altitude_km <= 2000:
            raise ValueError("drag altitude limits must satisfy 150 <= reentry < maximum <= 2000 km")


def atmospheric_density(altitude_km):
    """Nominal density; clamp below 150 km, extrapolate last scale height above 1000 km."""
    h = np.asarray(altitude_km, dtype=float)
    log_rho = np.interp(h, _HEIGHT, _LOG_DENSITY)
    slope = (_LOG_DENSITY[-1] - _LOG_DENSITY[-2]) / (_HEIGHT[-1] - _HEIGHT[-2])
    return np.exp(log_rho + np.maximum(h - _HEIGHT[-1], 0) * slope)


def area_to_mass(am, radius_m, mass_kg):
    """Use measured/modelled fragment A/M, otherwise spherical projected area / mass."""
    am, radius_m, mass_kg = np.broadcast_arrays(am, radius_m, mass_kg)
    if np.any(~np.isfinite(mass_kg) | (mass_kg <= 0) | ~np.isfinite(radius_m) | (radius_m < 0)):
        raise ValueError("drag requires positive mass and finite nonnegative radius")
    if np.any(~np.isnan(am) & (~np.isfinite(am) | (am < 0))):
        raise ValueError("invalid area-to-mass ratio")
    return np.where(np.isnan(am), np.pi * radius_m**2 / mass_kg, am)


def decay_rates(a_km, e, am, config: DragConfig):
    """Orbit-averaged da/dt [km/s], de/dt [1/s] for drag opposite inertial velocity.

    Energy loss gives da/dt = -B*rho*v*a*(v^2*a/mu); angular momentum gives
    de/dt = -B*rho*v*(e+cos(nu)), with B=Cd*A/m, v in m/s in the common factor.
    Average over time using dM=(1-e*cos(E))*dE, not uniform true anomaly.
    Chunk quadrature to bound temporary memory on the full catalogue.
    """
    a, e, am = np.broadcast_arrays(np.asarray(a_km, float), np.asarray(e, float), np.asarray(am, float))
    shape = a.shape
    a, e, am = a.ravel(), e.ravel(), am.ravel()
    da, de = np.empty_like(a), np.empty_like(e)
    for start in range(0, len(a), 2048):
        s = slice(start, start + 2048)
        aa, ee = a[s, None], e[s, None]
        w = 1 - ee * _COS_E
        r = aa * w
        v = np.sqrt(orbits.MU * (2 / r - 1 / aa))
        cos_nu = (_COS_E - ee) / w
        factor = config.cd * am[s, None] * config.density_scale * atmospheric_density(r - orbits.R_MEAN) * v * 1000
        weights = w * _WEIGHTS
        da[s] = -np.sum(factor * aa * (v*v * aa / orbits.MU) * weights, axis=1)
        de[s] = -np.sum(factor * (ee + cos_nu) * weights, axis=1)
    # Spherical atmosphere circularizes; suppress quadrature roundoff near e=0.
    de = np.where(e <= 1e-10, 0., np.minimum(de, 0.))
    return da.reshape(shape), de.reshape(shape)


# Packed segment: start, end, a,e,i,raan,argp,m, da,de, RAANdot,argpdot,Mdot.
START, END, A, E, I, RAAN, ARGP, M, DA, DE, DRAAN, DARGP, DM = range(13)
WIDTH = 13


def segment_elements(segments, t):
    """Evaluate a segment; after its end extrapolate the final orbit without further drag."""
    dt = np.asarray(t) - segments[:, START]
    elapsed = np.minimum(dt, segments[:, END] - segments[:, START])
    a = segments[:, A] + elapsed * segments[:, DA]
    e = np.clip(segments[:, E] + elapsed * segments[:, DE], 0, 1 - 1e-12)
    el = orbits.Elements(a, e, segments[:, I].copy(),
                         segments[:, RAAN] + elapsed * segments[:, DRAAN],
                         segments[:, ARGP] + elapsed * segments[:, DARGP],
                         segments[:, M] + elapsed * segments[:, DM])
    extra = dt - elapsed
    if np.any(extra):
        el = orbits.rebase(el, extra)
    return el


def maintained_segments(el, t):
    """Constant a/e, J2 propagation, no legacy ndot; one unbounded history segment."""
    t = np.broadcast_to(t, el.a.shape)
    zero = np.zeros_like(el.a)
    return np.column_stack((t, np.full_like(el.a, np.inf), el.a, el.e, el.i,
                            el.raan, el.argp, el.m, zero, zero, *orbits.secular_rates(el)))


def plan_segments(el, t, am, config: DragConfig):
    """Plan adaptive linear a/e segments with midpoint rates, independent of sampling ticks."""
    t = np.broadcast_to(t, el.a.shape)
    da, de = decay_rates(el.a, el.e, am, config)
    dt = np.minimum(config.max_step_s, config.max_axis_change_km / np.maximum(-da, 1e-300))
    dt = np.minimum(dt, 0.002 / np.maximum(-de, 1e-300))
    mid_a = el.a + da * dt / 2
    mid_e = np.maximum(el.e + de * dt / 2, 0)
    da, de = decay_rates(mid_a, mid_e, am, config)
    # Bound the actual midpoint step as well, then maintain e>=0.
    dt = np.minimum(dt, config.max_axis_change_km / np.maximum(-da, 1e-300))
    dt = np.minimum(dt, 0.002 / np.maximum(-de, 1e-300))
    de = np.maximum(de, -el.e / dt)
    threshold = orbits.R_MEAN + config.reentry_altitude_km
    def perigee_after(d):
        return (el.a + da * d) * (1 - el.e - de * d)
    hit = perigee_after(dt) <= threshold
    if hit.any():
        lo, hi = np.zeros_like(dt), dt.copy()
        for _ in range(40):
            mid = (lo + hi) / 2
            below = perigee_after(mid) <= threshold
            hi = np.where(below, mid, hi)
            lo = np.where(below, lo, mid)
        dt = np.where(hit, hi, dt)
    midpoint = orbits.Elements(el.a + da * dt / 2, el.e + de * dt / 2, el.i,
                              el.raan, el.argp, el.m)
    rates = orbits.secular_rates(midpoint)
    if np.any(t + dt <= t):
        raise ValueError("drag step too small for floating-point simulation time")
    return np.column_stack((t, t + dt, el.a, el.e, el.i, el.raan, el.argp, el.m,
                            da, de, *rates))
