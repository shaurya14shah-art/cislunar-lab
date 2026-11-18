"""Newtonian gravitational field quantities for an arbitrary set of point masses.

Sign conventions (physics convention):
    potential  Phi(x) = - sum_j G M_j / |x - x_j|     (J/kg, negative, zero far away)
    field      g(x)   = - grad Phi                     (m/s^2, points towards mass)
    tides      T_ij   = - d2 Phi / dx_i dx_j           (1/s^2)

The tidal tensor answers: "two nearby free particles separated by a small
vector d accelerate apart by  T d."  Along the line to a single mass
T = +2GM/r^3 (stretching); across it T = -GM/r^3 (squeezing).  Its trace is
zero in vacuum (Laplace's equation), which is why tides stretch one way and
squeeze the others - they never just compress or just expand.
"""
from __future__ import annotations
import numpy as np
from ..cr3bp.lagrange import lagrange_points
from ..constants import C


def potential(x, gm, pos):
    """Phi at points x (shape (..., 3))."""
    x = np.asarray(x, float)
    out = np.zeros(x.shape[:-1])
    for g, p in zip(gm, pos):
        if g:
            out -= g / np.linalg.norm(x - np.asarray(p, float), axis=-1)
    return out


def field(x, gm, pos):
    """Gravitational acceleration g at points x (shape (..., 3))."""
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    for g, p in zip(gm, pos):
        if g:
            d = np.asarray(p, float) - x
            r = np.linalg.norm(d, axis=-1, keepdims=True)
            out += g * d / r ** 3
    return out


def tidal_tensor(x, gm, pos):
    """3x3 tidal tensor at a single point x:  T = sum GM (3 n n^T - I) / r^3."""
    x = np.asarray(x, float)
    T = np.zeros((3, 3))
    for g, p in zip(gm, pos):
        if g:
            d = x - np.asarray(p, float)
            r = np.linalg.norm(d)
            n = d / r
            T += g * (3 * np.outer(n, n) - np.eye(3)) / r ** 3
    return T


def tidal_principal_axes(x, gm, pos):
    """Eigen-decomposition of the tidal tensor: (values, axes) sorted from the
    strongest stretch to the strongest squeeze.  axes[:, k] is the k-th direction."""
    w, V = np.linalg.eigh(tidal_tensor(x, gm, pos))
    order = np.argsort(w)[::-1]
    return w[order], V[:, order]


def potential_sheet(gm, pos, center, half_width, n=96, normal="z", offset=0.0, floor=None):
    """Potential sampled on a square grid in a coordinate plane through ``center``
    (the 'rubber sheet' picture).  Values are clipped at ``floor`` so the wells do
    not become infinitely deep spikes - the clip level is returned so the display
    can say so honestly.

    Returns dict(u, v, phi (n,n), floor).  For normal='z', u,v are x,y.
    """
    u = np.linspace(-half_width, half_width, n)
    U, V = np.meshgrid(u, u)
    c = np.asarray(center, float)
    W = np.full_like(U, offset)
    if normal == "z":
        pts = np.stack([c[0] + U, c[1] + V, c[2] + W], -1)
    elif normal == "y":
        pts = np.stack([c[0] + U, c[1] + W, c[2] + V], -1)
    else:
        pts = np.stack([c[0] + W, c[1] + U, c[2] + V], -1)
    phi = potential(pts, gm, pos)
    if floor is None:
        floor = float(np.percentile(phi, 1.0))
    return {"u": u, "v": u, "phi": np.maximum(phi, floor), "floor": floor}


def field_lines(gm, pos, seeds, step, n_steps=400, radius=None):
    """Trace gravitational field lines (curves everywhere tangent to g).

    Each line is followed "downhill" (along g) from its seed until it reaches a
    body or runs out of steps.  Uses 4th-order Runge-Kutta on dx/ds = g/|g|, so
    the step is a fixed *length*, which keeps lines smooth near and far alike.
    """
    pos = np.asarray(pos, float).reshape(-1, 3)
    rad = np.zeros(len(pos)) if radius is None else np.asarray(radius, float)

    def direction(x):
        g = field(x, gm, pos)
        return g / max(np.linalg.norm(g), 1e-300)

    lines = []
    for s in seeds:
        x = np.array(s, float)
        pts = [x.copy()]
        for _ in range(n_steps):
            k1 = direction(x)
            k2 = direction(x + 0.5 * step * k1)
            k3 = direction(x + 0.5 * step * k2)
            k4 = direction(x + step * k3)
            x = x + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            pts.append(x.copy())
            if np.any(np.linalg.norm(pos - x, axis=1) < np.maximum(rad, step)):
                break
        lines.append(np.array(pts))
    return lines


def effective_potential_pair(x, gm, pos, vel, pair):
    """Effective potential in the frame co-rotating with bodies pair=(i,j):

        Phi_eff = Phi - (1/2) |Omega x (x - x_bary)|^2

    (gravity of *all* bodies, plus the centrifugal potential of the pair's
    rotation).  Its stationary points are the pair's Lagrange points.
    """
    i, j = pair
    gm = np.asarray(gm, float)
    pos = np.asarray(pos, float)
    vel = np.asarray(vel, float)
    mt = gm[i] + gm[j]
    xb = (gm[i] * pos[i] + gm[j] * pos[j]) / mt
    d, dv = pos[j] - pos[i], vel[j] - vel[i]
    Om = np.cross(d, dv) / (d @ d)
    rel = np.asarray(x, float) - xb
    cen = np.cross(np.broadcast_to(Om, rel.shape), rel)
    return potential(x, gm, pos) - 0.5 * (cen ** 2).sum(-1)


def hill_radius(a, gm_small, gm_big, e=0.0):
    """Hill sphere radius  r_H = a (1-e) (m / 3M)^(1/3): inside it the small
    body's gravity dominates the tide of the big one, so it can keep satellites."""
    return a * (1 - e) * (gm_small / (3 * gm_big)) ** (1 / 3)


def roche_limits(gm_big, gm_small, radius_small):
    """Distances from the big body inside which the small one is torn apart.

    rigid  d = R_small (2 M/m)^(1/3)      self-gravity vs tide for a solid body
    fluid  d = 2.44 R_small (M/m)^(1/3)   a liquid body deforms first, so it breaks sooner
    (the 2.44 factor is Roche's classical result for a synchronously rotating fluid).
    Written with masses instead of densities: (M/m)^(1/3) R_small = R_big (rho_big/rho_small)^(1/3).
    """
    q = gm_big / gm_small
    return {"rigid": radius_small * (2 * q) ** (1 / 3), "fluid": 2.44 * radius_small * q ** (1 / 3)}


def pair_lagrange_points(gm, pos, vel, pair):
    """Lagrange points of bodies (i, j), placed in the *inertial* scene.

    Assumes the pair is on a (near-)circular orbit - the only case in which
    Lagrange points are true equilibria.  The returned 'eccentricity' lets the
    caller warn when that assumption is poor.
    """
    i, j = pair
    gm = np.asarray(gm, float); pos = np.asarray(pos, float); vel = np.asarray(vel, float)
    if gm[i] < gm[j]:
        i, j = j, i                                   # i is the heavier one
    mt = gm[i] + gm[j]
    mu = gm[j] / mt
    d, dv = pos[j] - pos[i], vel[j] - vel[i]
    L = np.linalg.norm(d)
    ex = d / L
    h = np.cross(d, dv)
    ez = h / np.linalg.norm(h)
    ey = np.cross(ez, ex)
    xb = (gm[i] * pos[i] + gm[j] * pos[j]) / mt
    evec = np.cross(dv, h) / mt - ex                  # eccentricity vector of the relative orbit
    pts = {k: xb + L * (v[0] * ex + v[1] * ey + v[2] * ez) for k, v in lagrange_points(mu).items()}
    return {"mu": float(mu), "separation": float(L), "eccentricity": float(np.linalg.norm(evec)),
            "points": pts}


def pair_report(bodies, pair, c=C):
    """Everything the interface shows when two bodies are selected together."""
    i, j = pair
    gm = np.array([b.gm for b in bodies]); pos = np.array([b.pos for b in bodies])
    vel = np.array([b.vel for b in bodies])
    big, small = (i, j) if gm[i] >= gm[j] else (j, i)
    lp = pair_lagrange_points(gm, pos, vel, pair)
    a = lp["separation"]
    roche = roche_limits(gm[big], gm[small], bodies[small].radius) if gm[small] > 0 else None
    return {
        "big": bodies[big].name, "small": bodies[small].name, "separation": a, "mu": lp["mu"],
        "eccentricity": lp["eccentricity"],
        "hill_radius": hill_radius(a, gm[small], gm[big]),
        "roche": roche,
        "inside_roche_fluid": bool(roche and a < roche["fluid"]),
        "inside_roche_rigid": bool(roche and a < roche["rigid"]),
        "contact": bool(a < bodies[big].radius + bodies[small].radius),
        "rs_big": 2 * gm[big] / c ** 2, "rs_small": 2 * gm[small] / c ** 2,
        "lagrange": {k: v.tolist() for k, v in lp["points"].items()},
        "circular_ok": bool(lp["eccentricity"] < 0.05),
    }
