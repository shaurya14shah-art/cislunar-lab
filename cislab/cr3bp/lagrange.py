"""The five Lagrange (libration) points and their linear stability.

A Lagrange point is a place where a body at rest in the rotating frame stays at
rest: grad(Omega) = 0.

* L1, L2, L3 lie on the x axis ("collinear").  Their positions solve a quintic
  equation that has no closed form, so we solve it numerically (Brent's method
  on dOmega/dx, which is monotonic between the singularities).
* L4, L5 form equilateral triangles with the primaries: x = 1/2 - mu, y = +-sqrt(3)/2.
  That result is exact for every mu.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from .core import grad_omega, jacobian, hessian_omega, omega


def _dOmega_dx(x, mu):
    return grad_omega(x, 0.0, 0.0, mu)[0]


def lagrange_points(mu):
    """Return {'L1': (x,y,z), ... 'L5': ...} in barycentric rotating coordinates."""
    if not 0 < mu <= 0.5:
        raise ValueError("mu must lie in (0, 0.5]; swap the bodies if mu > 0.5")
    tiny = 1e-12
    # Bracket each collinear root between the singularities at x=-mu and x=1-mu.
    x1 = brentq(_dOmega_dx, -mu + tiny ** 0.5, 1 - mu - 1e-9, args=(mu,), xtol=1e-15, rtol=1e-15)
    x2 = brentq(_dOmega_dx, 1 - mu + 1e-9, 2.0, args=(mu,), xtol=1e-15, rtol=1e-15)
    x3 = brentq(_dOmega_dx, -2.0, -mu - 1e-9, args=(mu,), xtol=1e-15, rtol=1e-15)
    s = np.sqrt(3) / 2
    return {"L1": np.array([x1, 0, 0.0]), "L2": np.array([x2, 0, 0.0]), "L3": np.array([x3, 0, 0.0]),
            "L4": np.array([0.5 - mu, s, 0.0]), "L5": np.array([0.5 - mu, -s, 0.0])}


def collinear_gamma(mu, point="L2"):
    """Distance from the *nearer primary* to L1/L2 (to the big primary for L3).
    Richardson's expansion uses this gamma as its length scale."""
    L = lagrange_points(mu)[point]
    if point in ("L1", "L2"):
        return abs(L[0] - (1 - mu))
    return abs(L[0] + mu)


def linear_frequencies(mu, point):
    """Characteristic rates of the linearised motion at a collinear point.

    Returns (lam, nu, omega_z):
      lam      real growth rate of the unstable direction  (distance ~ e^{lam t})
      nu       in-plane oscillation frequency (the Lyapunov-orbit frequency)
      omega_z  out-of-plane oscillation frequency

    Derivation.  On the x axis the Hessian of Omega is diagonal with
        Omega_xx = 1 + 2 c2,  Omega_yy = 1 - c2,  Omega_zz = -c2,
        c2 = (1-mu)/r1^3 + mu/r2^3.
    z decouples:  z'' = -c2 z  ->  omega_z = sqrt(c2).
    The in-plane characteristic equation is
        s^4 + (2 - c2) s^2 + (1 + 2c2)(1 - c2) = 0,
    whose roots are s = +-lam and s = +-i nu.
    """
    L = lagrange_points(mu)[point]
    H = hessian_omega(L[0], L[1], L[2], mu)
    c2 = -H[2, 2]
    disc = np.sqrt(9 * c2 * c2 - 8 * c2)
    lam = np.sqrt((c2 - 2 + disc) / 2)
    nu = np.sqrt((2 - c2 + disc) / 2)
    return float(lam), float(nu), float(np.sqrt(c2))


def classify_point(mu, point):
    """Eigenvalues of the 6x6 Jacobian at a Lagrange point plus a verdict.

    Verdict rules (linear theory):
      * any eigenvalue with positive real part  -> 'unstable (saddle x centre x centre)'
      * all eigenvalues purely imaginary        -> 'linearly stable (centre x centre x centre)'
    For L4/L5 linear stability holds only when mu < mu_Routh = 0.0385209...
    (Gascheau/Routh criterion  27 mu (1-mu) < 1).
    """
    L = lagrange_points(mu)[point]
    s = np.array([L[0], L[1], L[2], 0, 0, 0.0])
    eig = np.linalg.eigvals(jacobian(s, mu))
    tol = 1e-9
    unstable = bool(np.any(eig.real > tol))
    n_real = int(np.sum(np.abs(eig.imag) < tol))
    if unstable and n_real == 2:
        verdict = "unstable: saddle x centre x centre (one growing, one decaying, two oscillations)"
    elif unstable:
        verdict = "unstable: complex-unstable spiral (mu above the Routh limit)"
    else:
        verdict = "linearly stable: centre x centre x centre (three oscillations)"
    return {"point": point, "position": L.tolist(), "jacobi": float(2 * omega(*L, mu)),
            "eigenvalues": [[float(e.real), float(e.imag)] for e in eig],
            "unstable": unstable, "verdict": verdict}


MU_ROUTH = 0.5 * (1 - np.sqrt(23.0 / 27.0))   # 0.03852...: L4/L5 stable below this
