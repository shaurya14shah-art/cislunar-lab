"""Hill regions and zero-velocity surfaces.

Because C = 2*Omega - v^2 is conserved and v^2 can never be negative, a
spacecraft with Jacobi constant C can only ever be where

        2*Omega(x, y, z) >= C.

The boundary 2*Omega = C is the *zero-velocity surface* (ZVS): the spacecraft
would have to stop there.  The allowed side is the *Hill region*.  As C drops
(more energy) "necks" open first at L1, then L2, then L3, and finally the
forbidden region vanishes at L4/L5 (C = 3 exactly with our constant).
"""
from __future__ import annotations
import numpy as np
from .core import omega
from .lagrange import lagrange_points


def critical_jacobi(mu):
    """Jacobi constant of each Lagrange point: the energy at which its gate opens."""
    return {k: float(2 * omega(v[0], v[1], v[2], mu)) for k, v in lagrange_points(mu).items()}


def hill_slice(mu, C, plane="xy", offset=0.0, lim=(-1.5, 1.5, -1.5, 1.5), n=400):
    """2*Omega - C on a planar grid.  Positive = allowed, negative = forbidden.

    plane 'xy' at z=offset, 'xz' at y=offset, or 'yz' at x=offset.
    Returns (u, v, W) with W shaped (n, n) for direct use in contour plots.
    """
    u = np.linspace(lim[0], lim[1], n)
    v = np.linspace(lim[2], lim[3], n)
    U, V = np.meshgrid(u, v)
    O = np.full_like(U, offset)
    if plane == "xy":
        W = 2 * omega(U, V, O, mu) - C
    elif plane == "xz":
        W = 2 * omega(U, O, V, mu) - C
    elif plane == "yz":
        W = 2 * omega(O, U, V, mu) - C
    else:
        raise ValueError("plane must be 'xy', 'xz' or 'yz'")
    return u, v, W


def zvs_volume(mu, C, lim=(-1.5, 1.5, -1.5, 1.5, -1.0, 1.0), n=(72, 72, 48)):
    """2*Omega - C sampled on a 3D grid, for iso-surface extraction (marching
    cubes) in the front end.  The values are clipped so the 1/r spikes at the
    primaries do not swamp the iso-surface interpolation."""
    x = np.linspace(lim[0], lim[1], n[0])
    y = np.linspace(lim[2], lim[3], n[1])
    z = np.linspace(lim[4], lim[5], n[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    W = np.clip(2 * omega(X, Y, Z, mu) - C, -2.0, 2.0)
    return x, y, z, W


def regime(mu, C):
    """Plain-language description of which gateways are open at this energy."""
    c = critical_jacobi(mu)
    if C > c["L1"]:
        return "All gateways closed: motion is trapped around one primary (or far outside)."
    if C > c["L2"]:
        return "L1 neck open: transfer between the primaries is possible, escape is not."
    if C > c["L3"]:
        return "L1 and L2 necks open: the spacecraft can pass the secondary and escape."
    if C > c["L4"]:
        return "L1, L2 and L3 open: only tadpole-shaped islands around L4/L5 are forbidden."
    return "No forbidden region in the plane: every planar point is energetically reachable."
