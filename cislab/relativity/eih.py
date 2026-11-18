"""User-facing helpers for the 1PN (Einstein-Infeld-Hoffmann) N-body dynamics.

The compiled force law itself lives in ``cislab.nbody.kernels`` (with a term-by-
term explanation); this module adds the experiments built on it.
"""
from __future__ import annotations
import numpy as np
from ..constants import C, ARCSEC, JULIAN_CENTURY
from ..nbody.system import NBodySystem
from ..constants import Body


def newton_vs_1pn(system: NBodySystem, t_end, n_frames=600, rtol=1e-12):
    """Integrate the same initial conditions twice (Newtonian and 1PN) and
    return both histories plus the separation between them for every body."""
    hN = system.run(t_end, method="dop853", pn=False, n_frames=n_frames, rtol=rtol, update=False)
    hP = system.run(t_end, method="dop853", pn=True, n_frames=n_frames, rtol=rtol, update=False)
    div = np.linalg.norm(hP["pos"] - hN["pos"], axis=2)        # (frames, bodies), metres
    return {"newton": hN, "pn": hP, "divergence": div, "t": hN["t"]}


def perihelion_precession(gm_primary, gm_secondary, a, e, n_orbits=12, c=C, pn=True, rtol=1e-13):
    """Measure the perihelion advance of a two-body orbit by direct integration.

    Method: start at perihelion, integrate, and record the direction of the
    relative position vector at every later perihelion (detected as the instant
    when r.v changes sign from - to +).  A straight-line fit of that direction
    against orbit number gives the advance per orbit.
    Returns radians per orbit.
    """
    from scipy.integrate import solve_ivp
    from ..nbody.kernels import accel_newton, accel_eih_correction
    gm = np.array([gm_primary, gm_secondary], float)
    mt = gm.sum()
    rp = a * (1 - e)
    vp = np.sqrt(mt * (1 + e) / (a * (1 - e)))                 # vis-viva at perihelion
    f = gm[1] / mt
    pos0 = np.array([[-f * rp, 0, 0], [(1 - f) * rp, 0, 0]])
    vel0 = np.array([[0, -f * vp, 0], [0, (1 - f) * vp, 0]])
    period = 2 * np.pi * np.sqrt(a ** 3 / mt)

    def rhs(t, y):
        p, v = y[:6].reshape(2, 3), y[6:].reshape(2, 3)
        acc = np.empty((2, 3))
        accel_newton(p, gm, acc)
        if pn:
            da = np.empty((2, 3))
            accel_eih_correction(p, v, gm, c, acc, da)
            acc = acc + da
        return np.concatenate([v.ravel(), acc.ravel()])

    def peri(t, y):
        d, dv = y[3:6] - y[0:3], y[9:12] - y[6:9]
        return d @ dv
    peri.direction = 1
    atol = np.concatenate([np.full(6, 1e-14 * a), np.full(6, 1e-14 * vp)])
    sol = solve_ivp(rhs, (0, (n_orbits + 0.25) * period), np.concatenate([pos0.ravel(), vel0.ravel()]),
                    method="DOP853", rtol=rtol, atol=atol, events=peri)
    Y = sol.y_events[0]
    ang = np.unwrap(np.arctan2(Y[:, 4] - Y[:, 1], Y[:, 3] - Y[:, 0]))
    k = np.arange(1, len(ang) + 1)
    slope = np.polyfit(k, ang, 1)[0]
    return float(slope)


def mercury_precession_arcsec_per_century(pn=True):
    from ..constants import GM_SUN, GM_MERCURY, MERCURY_A, MERCURY_E, MERCURY_PERIOD_D, DAY
    per_orbit = perihelion_precession(GM_SUN, GM_MERCURY, MERCURY_A, MERCURY_E, pn=pn)
    orbits_per_century = JULIAN_CENTURY / (MERCURY_PERIOD_D * DAY)
    return per_orbit * orbits_per_century / ARCSEC
