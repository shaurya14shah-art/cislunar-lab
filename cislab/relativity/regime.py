"""Regime meter: *how relativistic is each body right now?*

Two small numbers decide everything:

    eps  = Phi/c^2 = sum_j G M_j / (r_ij c^2)     "how deep in a gravity well"
    beta2 = (v/c)^2                               "how fast"

For bound orbits the two are comparable (virial theorem: v^2 ~ GM/r).

* Newtonian gravity keeps terms of order 1.
* The first post-Newtonian approximation (1PN, the EIH equations) keeps terms
  of order eps.  What it throws away is of order eps^2 (2PN) - and from 2.5PN on
  the system radiates gravitational waves, which no conservative model captures.

So the *relative error of 1PN is about eps* (the neglected 2PN term divided by
the kept 1PN term), and the bands below are defined on that basis.  They are
conventions, not sharp physical boundaries; the meter always shows the raw
numbers so you can judge for yourself.
"""
from __future__ import annotations
import numpy as np
from ..constants import C

BANDS = [  # (upper limit of x = max(eps, beta^2), key, message)
    (1e-8, "newtonian", "Newtonian gravity is accurate to better than 1 part in 10⁸; 1PN is a tiny, trustworthy correction."),
    (1e-3, "1pn-good", "Weak field: 1PN is trustworthy (its own error is below 0.1% of the correction)."),
    (3e-2, "1pn-ok", "Moderately relativistic: 1PN captures the leading effect, but 2PN terms are a few % of it."),
    (1e-1, "1pn-marginal", "WARNING: strong field. 1PN is only a rough guide here (errors of order 10% of the correction)."),
    (np.inf, "breakdown", "WARNING: 1PN has broken down. Only a full general-relativistic treatment is meaningful."),
]


def lorentz_gamma(v, c=C):
    """gamma = 1/sqrt(1 - v^2/c^2).  NaN for v >= c (not allowed for matter)."""
    b2 = (np.asarray(v, float) / c) ** 2
    return np.where(b2 < 1, 1.0 / np.sqrt(np.clip(1 - b2, 1e-300, None)), np.nan)


def kinetic_energy_ratio(v, c=C):
    """(relativistic kinetic energy) / (Newtonian kinetic energy)
       = (gamma - 1) c^2 / (v^2 / 2).   Equals 1 + (3/4)(v/c)^2 + ... for slow motion."""
    v = np.asarray(v, float)
    b2 = (v / c) ** 2
    g = lorentz_gamma(v, c)
    # for tiny b2, (gamma-1) loses all digits; use the series instead
    series = 1 + 0.75 * b2 + 0.625 * b2 ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        exact = (g - 1) / (0.5 * b2)
    return np.where(b2 < 1e-6, series, exact)


def band(x):
    for lim, key, msg in BANDS:
        if x < lim:
            return key, msg
    return BANDS[-1][1], BANDS[-1][2]


def regime_report(gm, pos, vel, radius, names=None, c=C):
    """Per-body regime numbers.  Velocities are measured in the frame supplied
    (use the barycentric frame for a meaningful answer)."""
    gm, pos, vel, radius = map(lambda a: np.asarray(a, float), (gm, pos, vel, radius))
    n = len(gm)
    out = []
    for i in range(n):
        eps_ext = 0.0
        for j in range(n):
            if j != i and gm[j] > 0:
                eps_ext += gm[j] / (np.linalg.norm(pos[i] - pos[j]) * c * c)
        v = float(np.linalg.norm(vel[i]))
        b2 = (v / c) ** 2
        compact = gm[i] / (radius[i] * c * c) if radius[i] > 0 else np.inf
        x = max(eps_ext, b2)
        key, msg = band(x)
        out.append({
            "name": names[i] if names else f"body {i}",
            "eps_orbit": eps_ext,                # potential of the OTHER bodies at this body
            "beta2": b2,
            "compactness": compact,              # GM/(R c^2) at the body's own surface
            "rs_over_R": 2 * compact,            # Schwarzschild radius / radius (1 = black hole)
            "gamma": float(lorentz_gamma(v, c)),
            "ke_ratio": float(kinetic_energy_ratio(v, c)),
            "pn1_fractional_size": x,            # size of the 1PN correction relative to Newton
            "pn1_fractional_error": x * x,       # size of what 1PN neglects, relative to Newton
            "band": key, "message": msg,
            "surface_warning": bool(compact > 1e-2),
        })
    return out
