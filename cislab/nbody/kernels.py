"""Compiled force laws: Newtonian gravity and the 1PN (Einstein-Infeld-Hoffmann) correction.

Everything here works with gravitational parameters gm_j = G*M_j (m^3/s^2 in SI),
because GM is what astronomers actually measure well; G alone is known to only
about 5 significant figures.
"""
from __future__ import annotations
import numpy as np
from .._jit import njit


@njit(cache=True)
def accel_newton(pos, gm, out):
    """a_i = sum_j gm_j (r_j - r_i) / |r_j - r_i|^3   (Newton's inverse-square law)."""
    n = pos.shape[0]
    for i in range(n):
        out[i, 0] = out[i, 1] = out[i, 2] = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            r2 = dx * dx + dy * dy + dz * dz
            inv = 1.0 / (r2 * np.sqrt(r2))
            out[i, 0] += gm[j] * dx * inv; out[i, 1] += gm[j] * dy * inv; out[i, 2] += gm[j] * dz * inv
            out[j, 0] -= gm[i] * dx * inv; out[j, 1] -= gm[i] * dy * inv; out[j, 2] -= gm[i] * dz * inv


@njit(cache=True)
def accel_eih_correction(pos, vel, gm, c, aN, out):
    """First post-Newtonian correction to every body's acceleration.

    Einstein, Infeld & Hoffmann (1938); written here in the form used for
    planetary ephemerides (Moyer 2003; reproduced as Eq. 2 of Holman et al. 2023
    with beta = gamma = 1).  ``aN`` is the Newtonian acceleration, used wherever
    the formula needs another body's acceleration (the error made is O(1/c^4)).

    For body i, with r_ij = |r_i - r_j| and mu_j = G M_j:

      c^2 * da_i = sum_j mu_j (r_j - r_i)/r_ij^3 *
                   [ - 4 sum_{l!=i} mu_l/r_il      <- i sits deep in everyone's potential
                     -   sum_{k!=j} mu_k/r_jk      <- j's own potential energy gravitates
                     + v_i^2 + 2 v_j^2 - 4 v_i.v_j <- kinetic energy / motion of the source
                     - 3/2 ((r_i - r_j).v_j / r_ij)^2
                     + 1/2 (r_j - r_i).a_j ]       <- the source is accelerating
                 + sum_j mu_j/r_ij^3 [(r_i - r_j).(4 v_i - 3 v_j)] (v_i - v_j)
                                                    <- velocity-dependent ("gravitomagnetic-like") push
                 + 7/2 sum_j mu_j a_j / r_ij        <- again from the source's acceleration

    Every term is ~ (v/c)^2 or GM/(r c^2) times the Newtonian force.
    """
    n = pos.shape[0]
    pot = np.zeros(n)                      # pot[i] = sum_{l != i} mu_l / r_il
    for i in range(n):
        for j in range(n):
            if i != j:
                dx = pos[i, 0] - pos[j, 0]; dy = pos[i, 1] - pos[j, 1]; dz = pos[i, 2] - pos[j, 2]
                pot[i] += gm[j] / np.sqrt(dx * dx + dy * dy + dz * dz)
    c2 = c * c
    for i in range(n):
        ax = ay = az = 0.0
        vi2 = vel[i, 0] ** 2 + vel[i, 1] ** 2 + vel[i, 2] ** 2
        for j in range(n):
            if j == i or gm[j] == 0.0:
                continue
            dx = pos[j, 0] - pos[i, 0]; dy = pos[j, 1] - pos[i, 1]; dz = pos[j, 2] - pos[i, 2]
            r = np.sqrt(dx * dx + dy * dy + dz * dz)
            r3 = r * r * r
            vj2 = vel[j, 0] ** 2 + vel[j, 1] ** 2 + vel[j, 2] ** 2
            vivj = vel[i, 0] * vel[j, 0] + vel[i, 1] * vel[j, 1] + vel[i, 2] * vel[j, 2]
            rdv = -(dx * vel[j, 0] + dy * vel[j, 1] + dz * vel[j, 2]) / r     # (r_i-r_j).v_j / r
            rda = dx * aN[j, 0] + dy * aN[j, 1] + dz * aN[j, 2]               # (r_j-r_i).a_j
            br = -4.0 * pot[i] - pot[j] + vi2 + 2.0 * vj2 - 4.0 * vivj - 1.5 * rdv * rdv + 0.5 * rda
            f = gm[j] / r3 * br
            ax += f * dx; ay += f * dy; az += f * dz
            s = -(dx * (4 * vel[i, 0] - 3 * vel[j, 0]) + dy * (4 * vel[i, 1] - 3 * vel[j, 1])
                  + dz * (4 * vel[i, 2] - 3 * vel[j, 2]))
            g = gm[j] / r3 * s
            ax += g * (vel[i, 0] - vel[j, 0]); ay += g * (vel[i, 1] - vel[j, 1]); az += g * (vel[i, 2] - vel[j, 2])
            q = 3.5 * gm[j] / r
            ax += q * aN[j, 0]; ay += q * aN[j, 1]; az += q * aN[j, 2]
        out[i, 0] = ax / c2; out[i, 1] = ay / c2; out[i, 2] = az / c2


@njit(cache=True)
def yoshida_run(pos, vel, gm, dt, nsteps, stride, cc, dd, radius, P, V, T, coll):
    """Fixed-step 4th-order Yoshida (symplectic) integration, Newtonian forces.

    Stores every ``stride``-th step in P, V, T.  ``coll[i, j]`` receives the first
    time bodies i and j overlapped (distance < R_i + R_j), or stays -1.
    Returns the number of stored frames.
    """
    n = pos.shape[0]
    a = np.empty((n, 3))
    m = 0
    P[0] = pos; V[0] = vel; T[0] = 0.0
    m = 1
    t = 0.0
    for step in range(1, nsteps + 1):
        for s in range(4):
            for i in range(n):
                for k in range(3):
                    pos[i, k] += cc[s] * dt * vel[i, k]
            if s < 3:
                accel_newton(pos, gm, a)
                for i in range(n):
                    for k in range(3):
                        vel[i, k] += dd[s] * dt * a[i, k]
        t += dt
        for i in range(n):
            for j in range(i + 1, n):
                if coll[i, j] < 0.0:
                    dx = pos[j, 0] - pos[i, 0]; dy = pos[j, 1] - pos[i, 1]; dz = pos[j, 2] - pos[i, 2]
                    if dx * dx + dy * dy + dz * dz < (radius[i] + radius[j]) ** 2:
                        coll[i, j] = t
        if (step % stride == 0 or step == nsteps) and m < P.shape[0]:
            P[m] = pos; V[m] = vel; T[m] = t      # the final step is always stored
            m += 1
    return m
