"""Weak-field light bending by any number of bodies (for the N-body scene).

In a weak field light behaves as if space had a refractive index
n = 1 + 2U/c^2 (U = sum GM/r > 0).  A ray therefore curves towards mass:

    d(k_hat)/ds = (2/c^2) * grad_perp(U)      (the component of grad U across the ray)

The factor 2 is Einstein's famous doubling of the Newtonian value: half comes
from time curvature (which Newton "knows about"), half from space curvature.
Valid only where U/c^2 << 1.  Near a compact object use the exact
``schwarzschild.light_ray`` instead; the front end switches automatically and
marks rays that enter a region where this approximation is not trustworthy.
"""
from __future__ import annotations
import numpy as np
from ..constants import C


def trace_ray(x0, k0, gm, pos, length, n=800, c=C, radius=None):
    """March a light ray through the field of point masses.

    x0, k0: start point and direction; length: path length to trace (m).
    Returns dict(points (m,3), max_eps, blocked) where max_eps is the largest
    U/c^2 met along the ray (the validity indicator) and blocked says whether
    the ray hit a body's surface.
    """
    gm = np.asarray(gm, float)
    pos = np.asarray(pos, float).reshape(-1, 3)
    x = np.array(x0, float)
    k = np.array(k0, float)
    k /= np.linalg.norm(k)
    ds = length / n
    pts = [x.copy()]
    max_eps, blocked = 0.0, False

    def bend(xx, kk):
        d = pos - xx
        r = np.linalg.norm(d, axis=1)
        g = (gm / r ** 3)[:, None] * d                 # grad U summed below
        gU = g.sum(0)
        return (2.0 / (c * c)) * (gU - (gU @ kk) * kk), float((gm / r).sum() / (c * c)), r

    for _ in range(n):                                  # midpoint (RK2) steps
        dk1, e1, r = bend(x, k)
        if radius is not None and np.any(r < np.asarray(radius)):
            blocked = True
            break
        xm = x + 0.5 * ds * k
        km = k + 0.5 * ds * dk1
        km /= np.linalg.norm(km)
        dk2, e2, _ = bend(xm, km)
        x = x + ds * km
        k = k + ds * dk2
        k /= np.linalg.norm(k)
        max_eps = max(max_eps, e1, e2)
        pts.append(x.copy())
    return {"points": np.array(pts), "max_eps": max_eps, "blocked": blocked, "direction": k}


def deflection_angle(b, gm, c=C):
    """Weak-field bending angle 4GM/(c^2 b) for impact parameter b."""
    return 4 * gm / (c * c * b)
