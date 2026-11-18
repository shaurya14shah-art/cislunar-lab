"""Invariant manifolds of a periodic orbit ("tubes").

A halo orbit near L1/L2 is unstable.  Its monodromy matrix M (the STM over one
period) has one eigenvalue lambda_u > 1 and its partner lambda_s = 1/lambda_u < 1.

* Push the spacecraft a hair along the *unstable* eigenvector and it drifts away
  exponentially: the set of all such departures is the unstable manifold W^u.
* The *stable* manifold W^s is the set of trajectories that fall ONTO the orbit.
  We find it by stepping a hair along the stable eigenvector and integrating
  **backwards in time**.  Flying W^s forward costs no fuel: it is a free
  "conveyor belt" into the halo, which is what ballistic capture exploits.

Because every point of the orbit has its own eigenvector, each manifold is a
two-dimensional tube in space.  The eigenvector at a later point is obtained by
transporting the one at t=0 with the STM:  v(t) = Phi(t, 0) v(0).

The step size epsilon
---------------------
x_start = x_orbit +- eps * v, with v scaled so that its *position part has unit
length*: eps is then literally the initial displacement in length units.
Too large and the linear approximation is wrong; too small and the trajectory
loiters for a long time before leaving.  ``convergence_check`` repeats the
computation at eps/10 and measures how far the tube moved.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from . import fastprop as fp
from .periodic import PeriodicOrbit
from ..constants import R_MOON_KM, R_EARTH_KM, LUNIT_EM_KM


def eigen_directions(orbit: PeriodicOrbit):
    """(lambda_u, v_u, lambda_s, v_s) at the orbit's reference point.

    Raises ValueError for (near-)stable orbits, which have no real manifold pair.
    """
    w, V = np.linalg.eig(orbit.monodromy)
    real = np.abs(w.imag) < 1e-8 * np.maximum(1.0, np.abs(w))
    wr = np.where(real, w.real, np.nan)
    iu = int(np.nanargmax(np.abs(wr)))
    if not np.abs(wr[iu]) > 1.0 + 1e-6:
        raise ValueError("orbit is linearly stable: it has no stable/unstable manifolds")
    is_ = int(np.nanargmin(np.abs(wr)))
    vu, vs = V[:, iu].real, V[:, is_].real
    return float(wr[iu]), vu / np.linalg.norm(vu[:3]), float(wr[is_]), vs / np.linalg.norm(vs[:3])


@dataclass
class ManifoldPoint:
    phase: float               # fraction of the period at which we left the orbit
    sign: int                  # +1 / -1 : which side of the orbit we stepped to
    t: np.ndarray              # times (negative for the stable manifold)
    y: np.ndarray              # states, shape (m, 6)
    crossing: np.ndarray | None  # state where the trajectory first met the section
    t_cross: float | None
    impact: bool               # hit the secondary before reaching the section


def _first_section_crossing(T, Y, mu, value, axis=0):
    """Refine the first crossing of the plane  state[axis] = value."""
    g = Y[:, axis] - value
    idx = np.where(g[:-1] * g[1:] < 0)[0]
    if len(idx) == 0:
        return None, None, None
    k = int(idx[0])
    f = lambda tt: fp.integrate(Y[k], T[k], tt, mu)[axis] - value
    tc = brentq(f, T[k], T[k + 1], xtol=1e-15, rtol=1e-15)
    return tc, fp.integrate(Y[k], T[k], tc, mu), k


def manifold(orbit: PeriodicOrbit, kind="stable", n_phases=40, eps=1e-6, sign=+1, t_max=8.0,
             section_x=None, r_secondary=R_MOON_KM / LUNIT_EM_KM, hmax=0.01, phases=None,
             uniform=True):
    """Grow one branch of the stable or unstable manifold.

    uniform   True (default): step by eps * |p(phase)| / max|p|, where p is the
              Floquet-periodic eigenvector, so that every start has the same Floquet
              amplitude.  The starts are then spread *evenly over the trajectories
              of the tube* and changing eps only slides them along it.
              False: the textbook recipe "step eps at every phase".  It draws the
              same tube but samples its trajectories unevenly, and the unevenness
              shifts with eps - which made pooled Monte Carlo results depend on eps
              (found by the robustness test; see docs/DECISIONS.md, entry 21).
    kind      'stable' (integrated backwards) or 'unstable' (forwards)
    sign      +1 or -1: the two branches (towards / away from the secondary)
    section_x if given, each trajectory is cut at its first crossing of x = section_x
    Returns a list of ManifoldPoint, one per phase.
    """
    mu = orbit.mu
    lu, vu, ls, vs = eigen_directions(orbit)
    v0 = vs if kind == "stable" else vu
    tdir = -1.0 if kind == "stable" else 1.0
    lam = ls if kind == "stable" else lu
    if phases is None:
        phases = np.arange(n_phases) / n_phases
    # Floquet-periodic vector p(phase) = Phi(phase) v0 * lam^(-phase): its position
    # length says how big the eigen-displacement is at each point of the orbit
    # once the overall exponential decay/growth has been divided out.
    if uniform:
        grid = np.arange(200) / 200.0
        pmax = 0.0
        for g in grid:
            Pg = np.eye(6) if g == 0 else fp.integrate_stm(orbit.state0, 0.0, g * orbit.period, mu)[1]
            pmax = max(pmax, np.linalg.norm((Pg @ v0)[:3]) * abs(lam) ** (-g))
    out = []
    for ph in phases:
        if ph == 0:
            s, Phi = orbit.state0.copy(), np.eye(6)
        else:
            s, Phi = fp.integrate_stm(orbit.state0, 0.0, ph * orbit.period, mu)
        v = Phi @ v0
        scale = np.linalg.norm(v[:3]) * abs(lam) ** (-ph) / pmax if uniform else 1.0
        v /= np.linalg.norm(v[:3])            # unit position length: the step is a distance
        start = s + sign * eps * scale * v
        T, Y = fp.trajectory(start, 0.0, tdir * t_max, mu, hmax=hmax)
        r2 = np.sqrt((Y[:, 0] - 1 + mu) ** 2 + Y[:, 1] ** 2 + Y[:, 2] ** 2)
        hit = np.where(r2 < r_secondary)[0]
        tc = sc = None
        impact = False
        if section_x is not None:
            # ignore the loitering phase next to the orbit: only look once we are
            # well clear of it (more than 1e3*eps away would be arbitrary; instead
            # we simply take the first genuine sign change of x - section_x).
            tc, sc, k = _first_section_crossing(T, Y, mu, section_x)
            if tc is not None:
                if len(hit) and hit[0] <= k:
                    impact, tc, sc = True, None, None
                    T, Y = T[:hit[0] + 1], Y[:hit[0] + 1]
                else:
                    T, Y = np.append(T[:k + 1], tc), np.vstack([Y[:k + 1], sc])
        if section_x is None and len(hit):
            impact = True
            T, Y = T[:hit[0] + 1], Y[:hit[0] + 1]
        out.append(ManifoldPoint(float(ph), sign, T, Y, sc, tc, impact))
    return out


def branch_towards_secondary(orbit, kind="stable", eps=1e-6):
    """Return the sign (+1/-1) of the branch that heads for the secondary body.

    Decided by actually flying one trajectory of each sign and seeing which one
    gets closer to the secondary - simpler and safer than reasoning about
    eigenvector orientation, which is arbitrary.
    """
    best, best_r = +1, np.inf
    for sg in (+1, -1):
        mp = manifold(orbit, kind, phases=[0.0], eps=eps, sign=sg, t_max=6.0)[0]
        r2 = np.min(np.sqrt((mp.y[:, 0] - 1 + orbit.mu) ** 2 + mp.y[:, 1] ** 2 + mp.y[:, 2] ** 2))
        if r2 < best_r:
            best, best_r = sg, r2
    return best


def section_curve(points):
    """Stack the section crossings of a manifold into an array (n, 6) + phases."""
    ok = [p for p in points if p.crossing is not None]
    return np.array([p.phase for p in ok]), np.array([p.crossing for p in ok])


def convergence_check(orbit, kind="stable", eps=1e-6, sign=+1, section_x=None, n_phases=60):
    """Compare the section curve at eps and eps/10.

    The two runs start at slightly different points, so individual trajectories
    are not comparable; what must agree is the *curve* the tube draws on the
    section.  We report the largest distance from a point of the eps curve to
    the (densely interpolated) eps/10 curve, in the 4 section coordinates
    (y, z, vy, vz).  Small compared with the tube's own size => converged.
    """
    if section_x is None:
        section_x = 1 - orbit.mu
    _, A = section_curve(manifold(orbit, kind, n_phases, eps, sign, section_x=section_x))
    _, B = section_curve(manifold(orbit, kind, 4 * n_phases, eps / 10, sign, section_x=section_x))
    if len(A) < 3 or len(B) < 3:
        return {"ok": False, "reason": "too few section crossings"}
    cols = [1, 2, 4, 5]
    a, b = A[:, cols], B[:, cols]
    # distance from each eps-point to the eps/10 *polyline* (point-to-segment)
    p0, p1 = b[:-1], b[1:]
    seg = p1 - p0
    L2 = np.maximum((seg ** 2).sum(-1), 1e-300)
    d = np.empty(len(a))
    for i, q in enumerate(a):
        tpar = np.clip(((q - p0) * seg).sum(-1) / L2, 0.0, 1.0)
        d[i] = np.sqrt((((p0 + tpar[:, None] * seg) - q) ** 2).sum(-1)).min()
    size = float(np.linalg.norm(np.ptp(b, axis=0)))
    return {"ok": True, "max_distance": float(d.max()), "median_distance": float(np.median(d)),
            "tube_size": size, "relative": float(d.max() / size), "eps": eps,
            "n_eps": len(A), "n_eps10": len(B)}
