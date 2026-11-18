"""Poincare sections: turning a tangle of 3D curves into a readable 2D picture.

Instead of drawing a whole trajectory we only mark the points where it pierces
a chosen plane (the "section") in a chosen direction.  A periodic orbit becomes
a few dots, a regular (quasi-periodic) orbit a closed curve, a chaotic orbit a
scattered cloud, and a manifold tube a closed loop whose inside is the set of
trajectories that will transit through the gateway.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from . import fastprop as fp
from .core import omega


def crossings(state0, mu, t_max, axis=1, value=0.0, direction=+1, max_cross=200, hmax=0.01):
    """All crossings of state[axis] = value with the given direction of travel.

    direction +1: the coordinate is increasing as it crosses; -1 decreasing; 0 both.
    Returns (times, states[n,6]).  Each crossing is refined to ~1e-14 by root
    finding on the compiled integrator.
    """
    T, Y = fp.trajectory(state0, 0.0, t_max, mu, hmax=hmax)
    g = Y[:, axis] - value
    idx = np.where(g[:-1] * g[1:] < 0)[0]
    ts, ss = [], []
    for k in idx[:max_cross]:
        if direction and np.sign(g[k + 1] - g[k]) != np.sign(direction):
            continue
        f = lambda tt: fp.integrate(Y[k], T[k], tt, mu)[axis] - value
        tc = brentq(f, T[k], T[k + 1], xtol=1e-15, rtol=1e-15)
        ts.append(tc)
        ss.append(fp.integrate(Y[k], T[k], tc, mu))
    return np.array(ts), np.array(ss).reshape(-1, 6)


def planar_map(mu, C, x_values, t_max=60.0, vx0=0.0, prograde=True):
    """Planar Poincare map on the section y = 0, vy > 0 at fixed Jacobi constant C.

    For every starting x (with y = 0, vx = vx0) the speed vy is *forced* by the
    Jacobi constant:  vy = sqrt(2*Omega - C - vx0^2).  Points whose square root
    would be negative lie in the forbidden region and are skipped.
    Returns a list of (x0, array of (x, vx) crossings).
    """
    out = []
    for x0 in x_values:
        v2 = 2 * omega(x0, 0.0, 0.0, mu) - C - vx0 * vx0
        if v2 <= 0:
            continue
        s0 = np.array([x0, 0.0, 0.0, vx0, (1 if prograde else -1) * np.sqrt(v2), 0.0])
        _, S = crossings(s0, mu, t_max, axis=1, value=0.0, direction=+1)
        if len(S):
            out.append((float(x0), S[:, [0, 3]]))
    return out
