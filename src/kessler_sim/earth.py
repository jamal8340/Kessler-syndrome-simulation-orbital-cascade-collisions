"""Earth's rotation angle, needed to orient the Earth model under the satellites.

Positions are in the TEME frame (as SGP4); the matching Earth rotation angle is the
Greenwich mean sidereal time (IAU-82), the same function SGP4 uses internally.
UTC is used in place of UT1 (they differ by < 0.9 s, i.e. < 0.004 degrees).
"""

from __future__ import annotations

from sgp4.propagation import gstime


def earth_rotation_angle(jd_utc: float) -> float:
    """Greenwich mean sidereal angle [rad, 0..2pi] at the given Julian date."""
    return float(gstime(jd_utc))


def earth_rotation_at(sim, t: float) -> float:
    """Earth rotation angle at sim time ``t`` [s since the scenario epoch]."""
    return earth_rotation_angle(sim.epoch_jd + t / 86400.0)
