"""Periodic orbits: differential correction, monodromy, stability, continuation.

Symmetry trick (Howell 1984)
----------------------------
The CR3BP is unchanged by the mirror  (x, y, z, t) -> (x, -y, z, -t).  So if a
trajectory leaves the x-z plane perpendicularly (vx = vz = 0, y = 0) and later
hits that plane perpendicularly again, the mirror image closes the loop: we
have a periodic orbit with period twice the travel time.

So an orbit is described by four numbers  X = (x0, z0, vy0, tau)  where tau is
the half period, and it must satisfy three equations

        F(X) = [ y(tau), vx(tau), vz(tau) ] = 0.

Three equations, four unknowns -> a one-parameter *family* of orbits (the halo
family).  Newton's method needs dF/dX, and that is exactly what the State
Transition Matrix provides:

        dF/d(x0, z0, vy0) = Phi[rows y,vx,vz ; columns x,z,vy]
        dF/d(tau)         = [vy, ax, az] at the end point.

To pick ONE orbit we add a fourth equation (e.g. "z0 equals this value").
To walk along the family we add the pseudo-arclength equation
"(X - X_predicted) . tangent = 0", which keeps working when the family curve
turns around - natural-parameter continuation would fail there.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from .core import rhs, propagate, propagate_stm, jacobi, RTOL, ATOL
from .fastprop import integrate_stm, integrate
from .richardson import Richardson

_ROWS = [1, 3, 5]      # y, vx, vz must vanish at the half period
_COLS = [0, 2, 4]      # free initial components x0, z0, vy0


@dataclass
class PeriodicOrbit:
    mu: float
    state0: np.ndarray            # perpendicular crossing of the x-z plane
    period: float
    monodromy: np.ndarray = field(repr=False, default=None)
    residual: float = 0.0         # |F| after correction (how well the orbit closes)
    iterations: int = 0

    @property
    def jacobi(self):
        return float(jacobi(self.state0, self.mu))

    @property
    def multipliers(self):
        """Eigenvalues of the monodromy matrix (Floquet multipliers).
        They come in reciprocal pairs; two of them are always exactly 1."""
        return np.linalg.eigvals(self.monodromy)

    @property
    def stability_index(self):
        """k = (|lambda_max| + 1/|lambda_max|)/2.  k = 1 means marginally stable;
        k >> 1 means a displacement along the unstable direction is multiplied by
        about 2k every revolution.  (Same definition as NASA/JPL's catalogue;
        verified numerically in the tests.)"""
        lmax = np.max(np.abs(self.multipliers))
        return float(0.5 * (lmax + 1 / lmax))

    @property
    def unstable_multiplier(self):
        return float(np.max(np.abs(self.multipliers)))

    def sample(self, n=400, rtol=RTOL, atol=ATOL):
        """n states equally spaced in time around the orbit: returns (t, states[6,n])."""
        t = np.linspace(0, self.period, n, endpoint=False)
        y = np.empty((6, n))
        s = self.state0.copy()
        y[:, 0] = s
        for i in range(1, n):                      # march from sample to sample
            s = integrate(s, t[i - 1], t[i], self.mu, rtol=rtol, atol=atol)
            y[:, i] = s
        return t, y

    def amplitudes(self, n=720):
        """(Ax, Ay, Az): half-extents in x and y and max |z|, non-dimensional."""
        _, y = self.sample(n)
        return (float(np.ptp(y[0]) / 2), float(np.ptp(y[1]) / 2), float(np.max(np.abs(y[2]))))

    def to_dict(self):
        ax, ay, az = self.amplitudes(360)
        return {"mu": self.mu, "state0": self.state0.tolist(), "period": self.period,
                "jacobi": self.jacobi, "stability_index": self.stability_index,
                "unstable_multiplier": self.unstable_multiplier, "Ax": ax, "Ay": ay, "Az": az,
                "residual": self.residual}


def _F_and_DF(X, mu, rtol=RTOL, atol=ATOL):
    """Shooting residual F(X) and its 3x4 Jacobian."""
    x0, z0, vy0, tau = X
    s0 = np.array([x0, 0.0, z0, 0.0, vy0, 0.0])
    # Compiled DOP853 (same algorithm/coefficients as SciPy; see fastprop.py).
    sf, Phi = integrate_stm(s0, 0.0, tau, mu, rtol=rtol, atol=atol)
    F = sf[_ROWS]
    DF = np.empty((3, 4))
    DF[:, :3] = Phi[np.ix_(_ROWS, _COLS)]
    DF[:, 3] = rhs(tau, sf, mu)[_ROWS]
    return F, DF


def _finish(X, mu, res, it, rtol=RTOL, atol=ATOL):
    s0 = np.array([X[0], 0.0, X[1], 0.0, X[2], 0.0])
    _, M = integrate_stm(s0, 0.0, 2 * X[3], mu, rtol=rtol, atol=atol)
    return PeriodicOrbit(mu, s0, 2 * X[3], M, res, it)


def correct(X, mu, extra_row, extra_rhs_fn, tol=1e-12, max_iter=40, rtol=RTOL, atol=ATOL):
    """Newton iteration on the square system  [F(X); g(X)] = 0  where the scalar
    side condition is  g(X) = extra_row . X - extra_rhs_fn()  (linear in X).

    Returns (X, residual_norm, iterations).  Raises RuntimeError if it does not
    converge - a failed correction must never be silently accepted.
    """
    X = np.array(X, float)
    for it in range(1, max_iter + 1):
        F, DF = _F_and_DF(X, mu, rtol, atol)
        g = float(extra_row @ X - extra_rhs_fn())
        res = float(np.linalg.norm(F))
        if res < tol and abs(g) < tol:
            return X, res, it
        J = np.vstack([DF, extra_row])
        dX = np.linalg.solve(J, -np.append(F, g))
        # Damp very large steps: far from the solution the linear model is poor.
        scale = min(1.0, 0.05 / max(np.max(np.abs(dX)), 1e-300))
        X = X + scale * dX
    raise RuntimeError(f"differential correction did not converge (|F| = {res:.2e})")


def correct_fixed(X, mu, fix="z0", **kw):
    """Correct to the family member with a prescribed x0, z0 or period."""
    idx = {"x0": 0, "z0": 1, "vy0": 2, "tau": 3}[fix]
    row = np.zeros(4)
    row[idx] = 1.0
    target = X[idx]
    Xc, res, it = correct(X, mu, row, lambda: target, **kw)
    return _finish(Xc, mu, res, it, kw.get("rtol", RTOL), kw.get("atol", ATOL)), Xc


def halo_from_richardson(mu, az, point="L2", north=True, **kw):
    """Richardson seed -> differentially corrected halo orbit of amplitude ~az.

    The seed keeps its z0 fixed during correction, so the corrected orbit has
    the same z at the plane crossing as the analytic guess.
    """
    R = Richardson(mu, point)
    s = R.state(az, 0.0, north)
    X = np.array([s[0], s[2], s[4], R.period(az / R.gamma) / 2])
    return correct_fixed(X, mu, "z0", **kw)


def lyapunov_orbit(mu, ax, point="L2", **kw):
    """Planar Lyapunov orbit of in-plane amplitude ~ax (non-dimensional) from the
    linear solution x = -Ax cos(nu t), y = k Ax sin(nu t)."""
    R = Richardson(mu, point)
    x_l = 1 - mu + (R.gamma if point == "L2" else -R.gamma)
    X = np.array([x_l - ax, 0.0, R.k * ax * R.lam, np.pi / R.lam])
    return correct_fixed(X, mu, "x0", **kw)


def tangent(X, mu, prev=None, rtol=RTOL, atol=ATOL):
    """Unit tangent of the family curve at X (null vector of the 3x4 Jacobian),
    oriented to agree with the previous tangent so we keep walking forward."""
    _, DF = _F_and_DF(X, mu, rtol, atol)
    _, _, Vt = np.linalg.svd(DF)
    tvec = Vt[-1]
    if prev is not None and tvec @ prev < 0:
        tvec = -tvec
    return tvec


def continue_family(X0, mu, n_steps, ds=0.005, direction=None, stop=None, verbose=False,
                    ds_min=1e-5, ds_max=0.02, rtol=RTOL, atol=ATOL):
    """Pseudo-arclength continuation from a converged X0.

    direction: optional 4-vector; the first tangent is flipped to have a positive
               dot product with it (e.g. [0,1,0,0] to start towards larger z0).
    stop:      optional function(orbit) -> bool to end early (e.g. on lunar impact).
    Returns a list of PeriodicOrbit, including the starting orbit.
    """
    X = np.array(X0, float)
    tvec = tangent(X, mu, direction, rtol, atol)
    orbits = [_finish(X, mu, 0.0, 0, rtol, atol)]
    step = ds
    while len(orbits) <= n_steps:
        Xp = X + step * tvec
        try:
            Xn, res, it = correct(Xp, mu, tvec, lambda: float(tvec @ Xp), rtol=rtol, atol=atol,
                                  max_iter=12)
        except (RuntimeError, np.linalg.LinAlgError):
            step *= 0.5                       # predictor too bold: retry with half the step
            if step < ds_min:
                break
            continue
        tnew = tangent(Xn, mu, tvec, rtol, atol)
        X, tvec = Xn, tnew
        orb = _finish(X, mu, res, it, rtol, atol)
        orbits.append(orb)
        if verbose:
            print(f"  x0={X[0]:.6f} z0={X[1]:.6f} T={2*X[3]:.5f} k={orb.stability_index:.2f} it={it}")
        if stop is not None and stop(orb):
            break
        if it <= 3:
            step = min(step * 1.3, ds_max)    # easy step: be bolder next time
    return orbits
