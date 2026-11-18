"""Building and querying whole families of halo orbits.

Naming convention used throughout the project (same as NASA/JPL's catalogue):
an orbit is **northern** when z > 0 at the perpendicular x-z plane crossing that
lies on the *far* side of the libration point (away from the Moon for L2).
Richardson's own "class I" label refers to the near-side crossing, so the two
conventions are opposite for L2; ``build_family`` hides that detail.
"""
from __future__ import annotations
import numpy as np
from .fastprop import integrate
from .periodic import (PeriodicOrbit, halo_from_richardson, continue_family, correct_fixed)
from ..constants import R_MOON_KM, LUNIT_EM_KM


def far_side_X(orbit: PeriodicOrbit):
    """Describe the orbit by its *other* perpendicular crossing (half a period later)."""
    h = integrate(orbit.state0, 0.0, orbit.period / 2, orbit.mu)
    return np.array([h[0], h[2], h[4], orbit.period / 2])


def perilune_radius(orbit: PeriodicOrbit):
    _, y = orbit.sample(2000)
    return float(np.min(np.sqrt((y[0] - 1 + orbit.mu) ** 2 + y[1] ** 2 + y[2] ** 2)))


def build_family(mu, point="L2", north=True, seed_az=0.03, ds_max=0.01, max_steps=700,
                 min_period=1.2, r_body=R_MOON_KM / LUNIT_EM_KM):
    """Trace the halo family from near its birth (planar bifurcation) to the
    near-rectilinear end, stopping before the orbit would hit the secondary.

    Returns orbits ordered from small amplitude to large/NRHO, each stored at its
    far-side crossing so that state0 = (x_far, 0, z_far, 0, vy, 0).
    """
    # Richardson class with z<0 on the near side is "northern" in the JPL sense for L2.
    rich_north = (not north) if point == "L2" else north
    orb, X = halo_from_richardson(mu, seed_az, point, rich_north)
    sgn = np.sign(X[1])

    def too_close(o):
        return o.period < min_period or perilune_radius(o) < 1.02 * r_body

    up = continue_family(X, mu, max_steps, ds=0.004, direction=np.array([0, sgn, 0, 0.0]),
                         ds_max=ds_max, stop=too_close)
    down = continue_family(X, mu, 200, ds=0.002, direction=np.array([0, -sgn, 0, 0.0]),
                           ds_max=0.004, stop=lambda o: abs(o.state0[2]) < 2e-3)
    fam = down[::-1][:-1] + up
    out = []
    for o in fam:
        Xf = far_side_X(o)
        s0 = np.array([Xf[0], 0.0, Xf[1], 0.0, Xf[2], 0.0])
        # The monodromy matrix at a different point of the same orbit has the same
        # eigenvalues, so stability is unchanged; we recompute it for the new start.
        out.append(_restart(o, s0))
    return out


def _restart(o, s0):
    from .fastprop import integrate_stm
    _, M = integrate_stm(s0, 0.0, o.period, o.mu)
    return PeriodicOrbit(o.mu, s0, o.period, M, o.residual, o.iterations)


def member_with_z(family, z_target, branch="first"):
    """Family member whose far-side crossing has z = z_target.

    The far-side z is not monotonic along the whole family (it peaks and then
    falls on the way to the NRHOs), so ``branch`` selects 'first' (classic halo
    part, before the peak) or 'last' (near-rectilinear part).  A neighbouring
    member is used as the Newton guess and z is held fixed during correction,
    so the result is a genuine converged orbit, not an interpolation.
    """
    z = np.array([o.state0[2] for o in family])
    s = np.sign(z[np.argmax(np.abs(z))])
    zt = abs(z_target)
    az = np.abs(z)
    ipk = int(np.argmax(az))
    idx = range(0, ipk) if branch == "first" else range(len(family) - 2, ipk - 1, -1)
    for i in idx:
        lo, hi = sorted((az[i], az[i + 1]))
        if lo <= zt <= hi:
            a, b = family[i], family[i + 1]
            w = 0.0 if hi == lo else (zt - az[i]) / (az[i + 1] - az[i])
            Xa = np.array([a.state0[0], a.state0[2], a.state0[4], a.period / 2])
            Xb = np.array([b.state0[0], b.state0[2], b.state0[4], b.period / 2])
            X = (1 - w) * Xa + w * Xb
            X[1] = s * zt
            orb, _ = correct_fixed(X, family[0].mu, "z0")
            return orb
    raise ValueError(f"no family member with |z_far| = {zt} on the '{branch}' branch "
                     f"(range {az.min():.4f} .. {az.max():.4f})")


def family_table(family):
    """Compact numeric summary (one row per orbit) for plotting and the front end."""
    rows = []
    for o in family:
        rows.append([o.state0[0], o.state0[2], o.state0[4], o.period, o.jacobi, o.stability_index])
    return np.array(rows)
