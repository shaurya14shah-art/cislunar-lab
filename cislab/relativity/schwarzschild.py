"""Exact test-particle motion around one non-rotating mass (Schwarzschild geometry).

This is the opposite limit from 1PN N-body: *one* massive body, *arbitrarily
strong* field, and a companion so light it does not disturb the geometry.
Units inside this module are geometric: G = c = 1 and lengths are measured in
units of M (so the horizon is at r = 2, the photon sphere at 3, the ISCO at 6).
``to_si`` / ``from_si`` helpers convert.

Equatorial geodesics.  Two constants of motion exist: energy per unit mass E
and angular momentum per unit mass L.  The radial motion is that of a particle
in an effective potential

    (dr/dtau)^2 = E^2 - V(r),     V(r) = (1 - 2M/r)(kappa + L^2/r^2)

with kappa = 1 for matter and 0 for light.  Differentiating,

    d2r/dtau2 = -kappa M/r^2 + L^2/r^3 - 3 M L^2/r^4.

The first two terms are Newton (gravity + centrifugal barrier).  The last term
is new: an *attractive* 1/r^4 pull that beats the centrifugal barrier close in.
It causes perihelion precession, the innermost stable circular orbit, and the
bending of light.

    dphi/dtau = L/r^2,      dt/dtau = E/(1 - 2M/r)     (t = far-away observer's clock)
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq
from ..constants import C

HORIZON, PHOTON_SPHERE, ISCO = 2.0, 3.0, 6.0     # in units of GM/c^2


def length_unit(gm, c=C):
    """GM/c^2 in metres: the natural length of a mass (half its Schwarzschild radius)."""
    return gm / c ** 2


def time_unit(gm, c=C):
    return gm / c ** 3


def V_eff(r, L, kappa=1.0, M=1.0):
    return (1 - 2 * M / r) * (kappa + L * L / (r * r))


def circular_L2(r, M=1.0):
    """L^2 of a circular orbit of radius r (exists only for r > 3M)."""
    return M * r * r / (r - 3 * M)


def circular_E(r, M=1.0):
    return (1 - 2 * M / r) / np.sqrt(1 - 3 * M / r)


def isco_numeric(M=1.0):
    """Find the ISCO *numerically*, as the radius where circular orbits have the
    least angular momentum.  Below it, lowering r would require *more* L, so an
    inward nudge is not resisted: the orbit is unstable.  (Expected: 6M.)"""
    res = minimize_scalar(lambda r: circular_L2(r, M), bounds=(3.05 * M, 40 * M), method="bounded",
                          options={"xatol": 1e-12})
    return float(res.x)


def photon_sphere_numeric(M=1.0):
    """Radius of the circular light orbit = maximum of V for kappa = 0. (Expected: 3M.)"""
    res = minimize_scalar(lambda r: -V_eff(r, 1.0, 0.0, M), bounds=(2.05 * M, 20 * M), method="bounded",
                          options={"xatol": 1e-12})
    return float(res.x)


def _rhs(tau, y, L, kappa, E, M):
    t, r, phi, pr = y
    return [E / (1 - 2 * M / r), pr, L / (r * r),
            -kappa * M / (r * r) + L * L / r ** 3 - 3 * M * L * L / r ** 4]


def geodesic(r0, phi0, E, L, tau_max, kappa=1.0, inward=False, M=1.0, n=2000, rtol=1e-12):
    """Integrate an equatorial geodesic.  Returns dict with tau, t, r, phi, x, y
    and 'plunged' (True if it crossed r = 2M + small, where we stop: Schwarzschild
    coordinate time diverges at the horizon)."""
    pr2 = E * E - V_eff(r0, L, kappa, M)
    if pr2 < -1e-12:
        raise ValueError("these E, L cannot reach r0 (inside the forbidden region)")
    pr0 = (-1 if inward else 1) * np.sqrt(max(pr2, 0.0))
    horizon = lambda tau, y, *a: y[1] - 2.02 * M
    horizon.terminal = True
    sol = solve_ivp(_rhs, (0, tau_max), [0.0, r0, phi0, pr0], args=(L, kappa, E, M), method="DOP853",
                    rtol=rtol, atol=1e-14, t_eval=np.linspace(0, tau_max, n), events=horizon)
    t, r, phi, _ = sol.y
    return {"tau": sol.t, "t": t, "r": r, "phi": phi, "x": r * np.cos(phi), "y": r * np.sin(phi),
            "plunged": bool(sol.status == 1)}


def bound_orbit(r_peri, r_apo, n_orbits=3, M=1.0, n=4000):
    """Bound orbit between two turning points; returns the geodesic and its
    measured precession per orbit.

    E and L follow from requiring dr/dtau = 0 at both turning points:
    E^2 = V(r_peri) = V(r_apo), two equations for two unknowns.
    """
    a1, a2 = 1 - 2 * M / r_peri, 1 - 2 * M / r_apo
    L2 = (a2 - a1) / (a1 / r_peri ** 2 - a2 / r_apo ** 2)
    if L2 <= 0:
        raise ValueError("no bound orbit with these turning points")
    L = np.sqrt(L2)
    E = np.sqrt(V_eff(r_peri, L, 1.0, M))
    prec = precession_per_orbit(r_peri, r_apo, M)
    a = 0.5 * (r_peri + r_apo)
    tau_orb = 2 * np.pi * np.sqrt(a ** 3 / M)            # Kepler estimate of one period
    g = geodesic(r_peri * (1 + 1e-12), 0.0, E, L, n_orbits * tau_orb * 1.05, M=M, n=n)
    g.update({"E": float(E), "L": float(L), "precession_per_orbit": prec})
    return g


def precession_per_orbit(r_peri, r_apo, M=1.0):
    """Exact perihelion advance per radial period (radians), by quadrature.

    phi swept between two perihelia = 2 * integral_{r_p}^{r_a} (L/r^2) / sqrt(E^2 - V) dr.
    With u = 1/r the integrand becomes du / sqrt(2M (u - u1)(u - u2)(u - u3)); the
    substitution u = u1 + (u2-u1) sin^2(x) removes the square-root singularities
    at both turning points so ordinary quadrature is accurate.
    """
    from scipy.integrate import quad
    u1, u2 = 1 / r_apo, 1 / r_peri
    u3 = 1 / (2 * M) - u1 - u2          # third root of the cubic (sum of roots = 1/2M)
    f = lambda x: 2.0 / np.sqrt(2 * M * (u3 - u1 - (u2 - u1) * np.sin(x) ** 2))
    val, _ = quad(f, 0, np.pi / 2, epsabs=1e-13, epsrel=1e-13)
    return float(2 * val - 2 * np.pi)


def precession_weak_field(a, e, gm, c=C):
    """Leading-order formula 6 pi G M / (c^2 a (1 - e^2)) radians per orbit."""
    return 6 * np.pi * gm / (c * c * a * (1 - e * e))


def light_deflection(b, M=1.0):
    """Total bending angle (radians) of a light ray with impact parameter b.

    Uses the orbit equation for light  d2u/dphi2 + u = 3 M u^2  (u = 1/r).
    Starting from infinity (u = 0, du/dphi = 1/b) the ray reaches closest
    approach and returns to u = 0 after sweeping an angle pi + deflection.
    Returns NaN if the ray is captured (b below the critical value 3*sqrt(3) M).
    """
    if b <= 3 * np.sqrt(3) * M:
        return float("nan")
    # Work with w = b*u so the variable is of order one (good for error control):
    #   w'' + w = 3 (M/b) w^2,   w(0) = 0, w'(0) = 1.
    k = 3 * M / b
    rhs = lambda phi, y: [y[1], -y[0] + k * y[0] ** 2]
    back = lambda phi, y: y[0]
    back.terminal, back.direction = True, -1
    sol = solve_ivp(rhs, (0, 4 * np.pi), [0.0, 1.0], method="DOP853", rtol=1e-13, atol=1e-16,
                    events=back, max_step=0.05)
    if not len(sol.t_events[0]):
        return float("nan")
    return float(sol.t_events[0][0] - np.pi)


def light_deflection_weak(b, M=1.0):
    """Einstein's weak-field result 4GM/(c^2 b) (twice the Newtonian 'corpuscle' value)."""
    return 4 * M / b


def light_ray(x0, y0, direction_angle, M=1.0, lam_max=200.0, n=1500):
    """Trace a photon in the equatorial plane from (x0, y0) heading along the
    given angle (radians from +x), as seen by a far-away static observer's
    coordinates.  Returns x, y arrays and 'captured'.

    Note: we set the initial coordinate direction equal to the requested one.
    Close to the hole a *local* observer would measure a slightly different
    angle; for the visualiser's purpose (rays launched from far away) this
    distinction is negligible and is documented rather than hidden.
    """
    r0 = np.hypot(x0, y0)
    phi0 = np.arctan2(y0, x0)
    dx, dy = np.cos(direction_angle), np.sin(direction_angle)
    dr = dx * np.cos(phi0) + dy * np.sin(phi0)                 # radial component
    dphi = (-dx * np.sin(phi0) + dy * np.cos(phi0)) / r0      # angular rate per unit length
    L = r0 * r0 * dphi
    f = 1 - 2 * M / r0
    E = np.sqrt(dr * dr + f * L * L / (r0 * r0))
    horizon = lambda lam, y, *a: y[1] - 2.02 * M
    horizon.terminal = True
    far = lambda lam, y, *a: y[1] - max(3 * r0, 60 * M)
    far.terminal = True
    sol = solve_ivp(_rhs, (0, lam_max), [0.0, r0, phi0, dr], args=(L, 0.0, E, M), method="DOP853",
                    rtol=1e-11, atol=1e-13, events=[horizon, far], dense_output=True)
    lam = np.linspace(0, sol.t[-1], n)
    y = sol.sol(lam)
    return {"x": y[1] * np.cos(y[2]), "y": y[1] * np.sin(y[2]), "r": y[1],
            "captured": bool(len(sol.t_events[0]) > 0)}


def clock_rate_static(r, M=1.0):
    """d(tau)/dt for a clock hovering at radius r:  sqrt(1 - 2M/r)."""
    return np.sqrt(1 - 2 * M / np.asarray(r, float))


def clock_rate_circular(r, M=1.0):
    """d(tau)/dt for a clock on a circular orbit:  sqrt(1 - 3M/r).
    (2M/r from gravity, a further M/r from orbital speed.)"""
    return np.sqrt(1 - 3 * M / np.asarray(r, float))
