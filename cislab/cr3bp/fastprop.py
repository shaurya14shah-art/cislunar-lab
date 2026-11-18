"""A compiled (Numba) copy of the DOP853 integrator for the CR3BP.

Why does this exist when SciPy already has DOP853?
SciPy's stepper is written in Python, so every step costs ~100 microseconds of
interpreter overhead.  The Monte Carlo experiment needs hundreds of thousands of
trajectories and the family continuation needs thousands of 42-state
propagations; compiled code makes those ~100x faster.

It is the *same algorithm with the same coefficients* (they are imported from
SciPy), and tests/test_fastprop.py checks that the two agree to ~1e-12.  The
SciPy path (cislab.cr3bp.core.propagate) remains the reference implementation
and is what handles general event detection.

Optional 1PN term
-----------------
With ``pn=1`` the right-hand side adds the first post-Newtonian acceleration of
a test particle in the field of the two primaries (Einstein-Infeld-Hoffmann
equations in the test-mass limit, see cislab.relativity.eih for the full
N-body version and the meaning of each term).  The primaries are kept on their
Newtonian circular orbit: their own 1PN correction only rescales the relation
between separation and mean motion by ~1e-11, which is absorbed when the
*observed* distance and month are used as units.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate._ivp import dop853_coefficients as _c
from .._jit import njit

_A = np.ascontiguousarray(_c.A[:12, :12])
_B = np.ascontiguousarray(_c.B)
_C = np.ascontiguousarray(_c.C[:12])
_E3 = np.ascontiguousarray(_c.E3)
_E5 = np.ascontiguousarray(_c.E5)


@njit(cache=True)
def pn_test_particle(x, y, z, vx, vy, vz, mu, c):
    """1PN (EIH) correction to the acceleration of a massless particle, in the
    rotating frame at the instant when it coincides with the inertial frame.

    Inputs are rotating-frame position/velocity; G = 1; c is the speed of light
    in CR3BP units.  Returns (ax, ay, az) to be *added* to the Newtonian terms.
    """
    # inertial velocity of the particle: v_rot + omega x r, omega = (0,0,1)
    ux, uy, uz = vx - y, vy + x, vz
    v2 = ux * ux + uy * uy + uz * uz
    # primaries: positions, inertial velocities, accelerations on the circular orbit
    px = (-mu, 1.0 - mu)
    pvy = (-mu, 1.0 - mu)
    pax = (mu, -(1.0 - mu))
    gm = (1.0 - mu, mu)
    d = 1.0                                  # separation of the primaries
    r = (np.sqrt((x - px[0]) ** 2 + y * y + z * z), np.sqrt((x - px[1]) ** 2 + y * y + z * z))
    pot_i = gm[0] / r[0] + gm[1] / r[1]      # sum_l mu_l / r_il
    ax = ay = az = 0.0
    for j in range(2):
        o = 1 - j
        dx, dy, dz = px[j] - x, -y, -z       # r_j - r_i
        r3 = r[j] ** 3
        vj2 = pvy[j] * pvy[j]
        vivj = uy * pvy[j]
        # (r_i - r_j) . v_j / r_ij
        rdotvj = (-dx * 0.0 + (-dy) * pvy[j]) / r[j]
        # (r_j - r_i) . a_j
        rdotaj = dx * pax[j]
        bracket = (-4.0 * pot_i - gm[o] / d + v2 + 2.0 * vj2 - 4.0 * vivj
                   - 1.5 * rdotvj * rdotvj + 0.5 * rdotaj)
        f = gm[j] / r3 * bracket
        ax += f * dx
        ay += f * dy
        az += f * dz
        # velocity-dependent term: {(r_i - r_j).(4 v_i - 3 v_j)} (v_i - v_j)
        s = (-dx) * (4.0 * ux) + (-dy) * (4.0 * uy - 3.0 * pvy[j]) + (-dz) * (4.0 * uz)
        g = gm[j] / r3 * s
        ax += g * ux
        ay += g * (uy - pvy[j])
        az += g * uz
        # (7/2) mu_j a_j / r_ij
        ax += 3.5 * gm[j] * pax[j] / r[j]
    c2 = c * c
    return ax / c2, ay / c2, az / c2


@njit(cache=True)
def _rhs(y, mu, n, pn, c, out):
    X, Y, Z, vx, vy, vz = y[0], y[1], y[2], y[3], y[4], y[5]
    X1, X2 = X + mu, X - 1.0 + mu
    r1s = X1 * X1 + Y * Y + Z * Z
    r2s = X2 * X2 + Y * Y + Z * Z
    r1, r2 = np.sqrt(r1s), np.sqrt(r2s)
    a3, b3 = (1.0 - mu) / (r1s * r1), mu / (r2s * r2)
    out[0], out[1], out[2] = vx, vy, vz
    out[3] = 2.0 * vy + X - a3 * X1 - b3 * X2
    out[4] = -2.0 * vx + Y - (a3 + b3) * Y
    out[5] = -(a3 + b3) * Z
    if pn == 1:
        px, py, pz = pn_test_particle(X, Y, Z, vx, vy, vz, mu, c)
        out[3] += px
        out[4] += py
        out[5] += pz
    if n == 42:
        a5, b5 = 3.0 * (1.0 - mu) / (r1s * r1s * r1), 3.0 * mu / (r2s * r2s * r2)
        hxx = 1.0 - a3 - b3 + a5 * X1 * X1 + b5 * X2 * X2
        hyy = 1.0 - a3 - b3 + (a5 + b5) * Y * Y
        hzz = -a3 - b3 + (a5 + b5) * Z * Z
        hxy = (a5 * X1 + b5 * X2) * Y
        hxz = (a5 * X1 + b5 * X2) * Z
        hyz = (a5 + b5) * Y * Z
        for j in range(6):
            p0, p1, p2 = y[6 + j], y[12 + j], y[18 + j]
            p3, p4, p5 = y[24 + j], y[30 + j], y[36 + j]
            out[6 + j], out[12 + j], out[18 + j] = p3, p4, p5
            out[24 + j] = hxx * p0 + hxy * p1 + hxz * p2 + 2.0 * p4
            out[30 + j] = hxy * p0 + hyy * p1 + hyz * p2 - 2.0 * p3
            out[36 + j] = hxz * p0 + hyz * p1 + hzz * p2


@njit(cache=True)
def dop853_step(y, f, h, mu, n, pn, c, rtol, atol, K, ynew, fnew, A, B, C, E3, E5):
    """One trial step.  Returns the scaled error norm (accept if <= 1)."""
    K[0, :n] = f[:n]
    for s in range(1, 12):
        for i in range(n):
            acc = 0.0
            for q in range(s):
                acc += A[s, q] * K[q, i]
            ynew[i] = y[i] + h * acc
        _rhs(ynew, mu, n, pn, c, K[s])
    for i in range(n):
        acc = 0.0
        for q in range(12):
            acc += B[q] * K[q, i]
        ynew[i] = y[i] + h * acc
    _rhs(ynew, mu, n, pn, c, fnew)
    K[12, :n] = fnew[:n]
    e5 = 0.0
    e3 = 0.0
    for i in range(n):
        sc = atol + rtol * max(abs(y[i]), abs(ynew[i]))
        a5 = 0.0
        a3 = 0.0
        for q in range(13):
            a5 += E5[q] * K[q, i]
            a3 += E3[q] * K[q, i]
        e5 += (a5 / sc) ** 2
        e3 += (a3 / sc) ** 2
    if e5 == 0.0 and e3 == 0.0:
        return 0.0
    return abs(h) * e5 / np.sqrt((e5 + 0.01 * e3) * n)


@njit(cache=True)
def _integrate(y0, t0, t1, mu, n, pn, c, rtol, atol, hmax, A, B, C, E3, E5):
    y = y0.copy()
    f = np.empty(n)
    fnew = np.empty(n)
    ynew = np.empty(n)
    K = np.empty((13, n))
    _rhs(y, mu, n, pn, c, f)
    direction = 1.0 if t1 >= t0 else -1.0
    h = direction * min(1e-3, hmax, abs(t1 - t0))
    t = t0
    nsteps = 0
    while (t1 - t) * direction > 0.0:
        if (t + h - t1) * direction > 0.0:
            h = t1 - t
        err = dop853_step(y, f, h, mu, n, pn, c, rtol, atol, K, ynew, fnew, A, B, C, E3, E5)
        if err <= 1.0:
            t += h
            for i in range(n):
                y[i] = ynew[i]
                f[i] = fnew[i]
            nsteps += 1
            fac = 10.0 if err == 0.0 else min(10.0, 0.9 * err ** (-0.125))
            h = direction * min(abs(h) * fac, hmax)
        else:
            h *= max(0.2, 0.9 * err ** (-0.125))
        if abs(h) < 1e-14:
            break                      # step-size underflow: we hit a singularity
    return y, nsteps


def integrate(y0, t0, t1, mu, rtol=3e-14, atol=1e-16, pn=False, c=np.inf, hmax=np.inf):
    """Propagate a 6-state from t0 to t1 (either direction).  Returns final state."""
    y, _ = _integrate(np.ascontiguousarray(y0, dtype=np.float64), float(t0), float(t1), float(mu), 6,
                      1 if pn else 0, float(c), rtol, atol, hmax, _A, _B, _C, _E3, _E5)
    return y


def integrate_stm(y0, t0, t1, mu, rtol=3e-14, atol=1e-16):
    """Propagate state + STM.  Returns (final 6-state, 6x6 Phi)."""
    y = np.concatenate([np.asarray(y0, float), np.eye(6).ravel()])
    y, _ = _integrate(y, float(t0), float(t1), float(mu), 42, 0, np.inf, rtol, atol, np.inf,
                      _A, _B, _C, _E3, _E5)
    return y[:6], y[6:].reshape(6, 6).copy()


@njit(cache=True)
def _trajectory(y0, t0, t1, mu, pn, c, rtol, atol, hmax, maxpts, A, B, C, E3, E5):
    n = 6
    y = y0.copy()
    f = np.empty(n); fnew = np.empty(n); ynew = np.empty(n)
    K = np.empty((13, n))
    T = np.empty(maxpts)
    Y = np.empty((maxpts, n))
    _rhs(y, mu, n, pn, c, f)
    direction = 1.0 if t1 >= t0 else -1.0
    h = direction * min(1e-3, hmax, abs(t1 - t0))
    t = t0
    T[0] = t
    Y[0] = y
    m = 1
    while (t1 - t) * direction > 0.0 and m < maxpts:
        if (t + h - t1) * direction > 0.0:
            h = t1 - t
        err = dop853_step(y, f, h, mu, n, pn, c, rtol, atol, K, ynew, fnew, A, B, C, E3, E5)
        if err <= 1.0:
            t += h
            for i in range(n):
                y[i] = ynew[i]
                f[i] = fnew[i]
            T[m] = t
            Y[m] = y
            m += 1
            fac = 10.0 if err == 0.0 else min(10.0, 0.9 * err ** (-0.125))
            h = direction * min(abs(h) * fac, hmax)
        else:
            h *= max(0.2, 0.9 * err ** (-0.125))
        if abs(h) < 1e-14:
            break
    return T[:m], Y[:m]


def trajectory(y0, t0, t1, mu, rtol=3e-14, atol=1e-16, pn=False, c=np.inf, hmax=0.01,
               maxpts=200_000):
    """Propagate and keep every accepted step (capped at ``hmax`` so the output
    is smooth enough to draw).  Returns (t[m], states[m,6])."""
    return _trajectory(np.ascontiguousarray(y0, dtype=np.float64), float(t0), float(t1), float(mu),
                       1 if pn else 0, float(c), rtol, atol, hmax, maxpts, _A, _B, _C, _E3, _E5)
