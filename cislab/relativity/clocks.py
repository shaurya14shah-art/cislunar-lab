"""Clock rates, time dilation and gravitational redshift in weak fields.

To first post-Newtonian order a clock at position x moving with speed v ticks at

    d(tau)/dt = 1 - U(x)/c^2 - v^2/(2 c^2),        U = sum_j G M_j / |x - x_j|  (> 0)

relative to coordinate time t (the time of a far-away observer at rest).
* -U/c^2      : gravitational time dilation - clocks deeper in a well run slow;
* -v^2/(2c^2) : special-relativistic time dilation - moving clocks run slow.
Valid while both terms are << 1 (see the regime meter).

Gravitational redshift between two static clocks: light emitted at A and
received at B is shifted by  z = (rate_B / rate_A) - 1, so z > 0 (red) when
the light climbs out of a well.
"""
from __future__ import annotations
import numpy as np
from ..constants import (C, GM_EARTH, GPS_A, DAY)


def potential_over_c2(x, gm, pos, c=C, exclude=None):
    x = np.asarray(x, float)
    U = 0.0
    for j, (g, p) in enumerate(zip(gm, pos)):
        if j == exclude or g == 0:
            continue
        U += g / np.linalg.norm(x - np.asarray(p, float))
    return U / (c * c)


def clock_rate(x, v, gm, pos, c=C, exclude=None):
    """d(tau)/dt - 1  (a small negative number; multiply by 86400e6 for microseconds/day)."""
    v2 = float(np.dot(v, v))
    return -potential_over_c2(x, gm, pos, c, exclude) - v2 / (2 * c * c)


def body_clock_rates(gm, pos, vel, radius, c=C):
    """Rate offset of a clock riding on each body's *surface* (own potential
    GM/R included) as it moves with the body.  Returns d(tau)/dt - 1 per body.
    For compact objects the exact static factor sqrt(1 - 2GM/(R c^2)) is used
    for the body's own field instead of the weak-field term."""
    out = []
    for i in range(len(gm)):
        ext = clock_rate(pos[i], vel[i], gm, pos, c, exclude=i)
        own = 2 * gm[i] / (radius[i] * c * c) if radius[i] > 0 else 1.0
        own_rate = np.sqrt(max(1 - own, 0.0)) - 1.0          # exact for a static non-rotating body
        out.append(ext + own_rate)
    return np.array(out)


def redshift(rate_emit, rate_recv):
    """z = lambda_received/lambda_emitted - 1 from the two clock-rate offsets."""
    return (1 + rate_recv) / (1 + rate_emit) - 1


def gps_clock_offset(a=GPS_A, gm=GM_EARTH, c=C, a1=6.3781370e6, J2=1.0826300e-3,
                     omega_e=7.2921151467e-5):
    """Net rate of a GPS satellite clock relative to a clock on Earth's geoid.

    Satellite (circular orbit, radius a):  potential GM/a plus v^2/2 = GM/(2a)
        -> runs slow by (3/2) GM/(a c^2) relative to coordinate time.
    Ground clock on the rotating geoid: potential including Earth's oblateness
        (J2) and the rotation speed at the equator,
        Phi0/c^2 = GM/(a1 c^2) (1 + J2/2) + (omega_e a1)^2 / (2 c^2).
        The geoid is an equipotential of exactly this combined potential, so one
        value serves every sea-level clock.
    Net fractional offset = Phi0/c^2 - (3/2) GM/(a c^2)  (positive: satellite fast).

    Formula and constants follow N. Ashby, Living Rev. Relativity 6 (2003) 1,
    Eq. (35), whose value is 4.4647e-10.  The Earth constants a1, J2, omega_e are
    the WGS-84 values Ashby uses; they are quoted, not re-derived.
    """
    sat = 1.5 * gm / (a * c * c)
    ground = gm / (a1 * c * c) * (1 + J2 / 2) + (omega_e * a1) ** 2 / (2 * c * c)
    net = ground - sat
    return {"satellite_term": sat, "ground_term": ground, "fractional_offset": net,
            "microseconds_per_day": net * DAY * 1e6,
            "gravitational_only_us_per_day": (gm / (a1 * c * c) - gm / (a * c * c)) * DAY * 1e6,
            "velocity_only_us_per_day": -(gm / (2 * a * c * c) - (omega_e * a1) ** 2 / (2 * c * c)) * DAY * 1e6}
