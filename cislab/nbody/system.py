"""The N-body sandbox: an isolated set of gravitating bodies you can edit freely.

Integrators
-----------
* ``yoshida4``  fixed-step, 4th-order, *symplectic*.  Symplectic methods do not
  conserve energy exactly, but their energy error oscillates and never drifts,
  so they are the right tool for long runs (thousands of orbits).  They need the
  force to depend on positions only, so they are used for Newtonian gravity.
* ``dop853``    adaptive 8th-order Runge-Kutta.  Very accurate per step, handles
  close encounters by shrinking the step, and accepts velocity-dependent forces
  - so this is what the 1PN (relativistic) option uses.  Energy drifts slowly
  but monotonically; we always report it.

Collisions are *flagged, never merged*: a real merger needs hydrodynamics that
this code does not pretend to have.
"""
from __future__ import annotations
import copy
import numpy as np
from scipy.integrate import solve_ivp
from ..constants import Body, PRESETS, C, G
from .kernels import accel_newton, accel_eih_correction, yoshida_run


def yoshida_coefficients():
    """Drift (c) and kick (d) weights of Yoshida's 4th-order composition of
    three leapfrog steps.  Note the middle step runs *backwards* (w0 < 0): that
    is what cancels the 3rd-order error of the plain leapfrog."""
    w1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
    w0 = -(2.0 ** (1.0 / 3.0)) * w1
    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1 / 2])
    d = np.array([w1, w0, w1, 0.0])
    return c, d


class NBodySystem:
    def __init__(self, bodies=None, c=C):
        self.bodies: list[Body] = [copy.deepcopy(b) for b in (bodies or [])]
        self.c = c
        self.t = 0.0

    # ---- editing -----------------------------------------------------------
    def add(self, body: Body | str, pos=None, vel=None, name=None):
        b = copy.deepcopy(PRESETS[body]) if isinstance(body, str) else copy.deepcopy(body)
        if pos is not None:
            b.pos = tuple(float(x) for x in pos)
        if vel is not None:
            b.vel = tuple(float(x) for x in vel)
        if name:
            b.name = name
        self.bodies.append(b)
        return len(self.bodies) - 1

    def remove(self, index):
        return self.bodies.pop(index)

    def edit(self, index, **changes):
        for k, v in changes.items():
            if k == "mass":
                k, v = "gm", v * G
            if not hasattr(self.bodies[index], k):
                raise AttributeError(f"Body has no field '{k}'")
            setattr(self.bodies[index], k, tuple(v) if isinstance(v, (list, np.ndarray)) else v)

    def add_circular_satellite(self, body, primary_index, distance, direction=(1, 0, 0),
                               normal=(0, 0, 1), name=None):
        """Place ``body`` on a circular two-body orbit around an existing body.
        Speed from  v^2 = G(M + m)/r ; the primary is given the opposite momentum
        so the pair's centre of mass keeps its velocity."""
        b = copy.deepcopy(PRESETS[body]) if isinstance(body, str) else copy.deepcopy(body)
        P = self.bodies[primary_index]
        d = np.asarray(direction, float); d /= np.linalg.norm(d)
        nrm = np.asarray(normal, float); nrm /= np.linalg.norm(nrm)
        tang = np.cross(nrm, d)
        v_rel = np.sqrt((P.gm + b.gm) / distance) * tang
        f = b.gm / (P.gm + b.gm)
        ppos, pvel = np.array(P.pos), np.array(P.vel)
        b.pos = tuple(ppos + (1 - f) * distance * d)
        b.vel = tuple(pvel + (1 - f) * v_rel)
        P.pos = tuple(ppos - f * distance * d)
        P.vel = tuple(pvel - f * v_rel)
        if name:
            b.name = name
        self.bodies.append(b)
        return len(self.bodies) - 1

    # ---- array views -------------------------------------------------------
    @property
    def gm(self):
        return np.array([b.gm for b in self.bodies], dtype=float)

    @property
    def radius(self):
        return np.array([b.radius for b in self.bodies], dtype=float)

    @property
    def pos(self):
        return np.array([b.pos for b in self.bodies], dtype=float).reshape(-1, 3)

    @property
    def vel(self):
        return np.array([b.vel for b in self.bodies], dtype=float).reshape(-1, 3)

    def set_state(self, pos, vel):
        for b, p, v in zip(self.bodies, pos, vel):
            b.pos, b.vel = tuple(map(float, p)), tuple(map(float, v))

    # ---- physics -----------------------------------------------------------
    def acceleration(self, pos=None, vel=None, pn=False):
        pos = self.pos if pos is None else np.ascontiguousarray(pos, float)
        vel = self.vel if vel is None else np.ascontiguousarray(vel, float)
        gm = self.gm
        a = np.empty_like(pos)
        accel_newton(pos, gm, a)
        if pn:
            da = np.empty_like(pos)
            accel_eih_correction(pos, vel, gm, self.c, a, da)
            return a + da
        return a

    def energy(self, pos=None, vel=None, pn=False):
        """Total energy divided by G (so that it is built from GM values).

        Newtonian:  E = sum 1/2 m v^2 - sum_{pairs} G m_i m_j / r_ij.
        With pn=True the 1PN-conserved energy of the EIH equations is returned
        (Landau & Lifshitz, Classical Theory of Fields, section 106).  That
        formula was not read from the book while building this project; it is
        validated numerically instead - see tests/test_relativity.py, which
        checks that it is conserved ~(v/c)^2 times better than the Newtonian one.
        """
        pos = self.pos if pos is None else np.asarray(pos, float)
        vel = self.vel if vel is None else np.asarray(vel, float)
        gm = self.gm
        n = len(gm)
        v2 = (vel ** 2).sum(1)
        E = 0.5 * (gm * v2).sum()
        c2 = self.c ** 2
        rinv = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    rinv[i, j] = 1.0 / np.linalg.norm(pos[i] - pos[j])
        E -= 0.5 * (gm[:, None] * gm[None, :] * rinv).sum()
        if pn:
            E += 0.375 * (gm * v2 * v2).sum() / c2
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    nij = (pos[i] - pos[j]) * rinv[i, j]
                    E += 1.5 * gm[i] * gm[j] * v2[i] * rinv[i, j] / c2
                    E -= 0.25 * gm[i] * gm[j] * rinv[i, j] / c2 * (
                        7 * vel[i] @ vel[j] + (vel[i] @ nij) * (vel[j] @ nij))
                    for k in range(n):
                        if k != i:
                            E += 0.5 * gm[i] * gm[j] * gm[k] * rinv[i, j] * rinv[i, k] / c2
        return float(E)

    def momentum(self, vel=None):
        vel = self.vel if vel is None else vel
        return (self.gm[:, None] * vel).sum(0)

    def angular_momentum(self, pos=None, vel=None):
        pos = self.pos if pos is None else pos
        vel = self.vel if vel is None else vel
        return (self.gm[:, None] * np.cross(pos, vel)).sum(0)

    def barycenter(self, pos=None, vel=None):
        pos = self.pos if pos is None else pos
        vel = self.vel if vel is None else vel
        gm = self.gm
        return (gm[:, None] * pos).sum(0) / gm.sum(), (gm[:, None] * vel).sum(0) / gm.sum()

    def to_barycentric(self):
        """Shift so the centre of mass is at rest at the origin (an isolated
        system's natural frame: total momentum is then zero for all time)."""
        r0, v0 = self.barycenter()
        self.set_state(self.pos - r0, self.vel - v0)

    # ---- integration -------------------------------------------------------
    def run(self, t_end, dt=None, method="yoshida4", pn=False, n_frames=600, rtol=1e-12, atol=None,
            update=True):
        """Advance the system by ``t_end`` seconds and return a history dict.

        history keys: t (F,), pos (F,N,3), vel (F,N,3), energy_drift (F,) relative,
        collisions [ {i, j, t, names} ], method, pn.
        """
        n = len(self.bodies)
        if n == 0:
            raise ValueError("the system has no bodies")
        pos0, vel0, gm, rad = self.pos.copy(), self.vel.copy(), self.gm, self.radius
        if pn and method == "yoshida4":
            raise ValueError("1PN forces depend on velocity, which breaks the symplectic splitting; "
                             "use method='dop853' for relativistic runs")
        if method == "yoshida4":
            if dt is None:
                raise ValueError("yoshida4 needs a fixed time step dt")
            nsteps = max(1, int(round(abs(t_end) / dt)))
            stride = max(1, nsteps // n_frames)
            F = nsteps // stride + 2
            P, V, T = np.empty((F, n, 3)), np.empty((F, n, 3)), np.empty(F)
            coll = -np.ones((n, n))
            cc, dd = yoshida_coefficients()
            m = yoshida_run(pos0.copy(), vel0.copy(), gm, np.sign(t_end) * dt, nsteps, stride, cc, dd,
                            rad, P, V, T, coll)
            P, V, T = P[:m], V[:m], T[:m]
        elif method == "dop853":
            c = self.c

            def rhs(t, y):
                p = y[:3 * n].reshape(n, 3)
                v = y[3 * n:].reshape(n, 3)
                a = np.empty((n, 3))
                accel_newton(p, gm, a)
                if pn:
                    da = np.empty((n, 3))
                    accel_eih_correction(p, v, gm, c, a, da)
                    a = a + da
                return np.concatenate([v.ravel(), a.ravel()])

            y0 = np.concatenate([pos0.ravel(), vel0.ravel()])
            if atol is None:   # absolute tolerance scaled to the problem's size
                L = max(np.abs(pos0).max(), 1.0)
                Vs = max(np.abs(vel0).max(), 1e-30)
                atol = np.concatenate([np.full(3 * n, 1e-13 * L), np.full(3 * n, 1e-13 * Vs)])
            T = np.linspace(0.0, t_end, n_frames + 1)
            sol = solve_ivp(rhs, (0.0, t_end), y0, method="DOP853", t_eval=T, rtol=rtol, atol=atol)
            if not sol.success:
                raise RuntimeError("integration failed: " + sol.message)
            T = sol.t
            P = sol.y[:3 * n].T.reshape(-1, n, 3)
            V = sol.y[3 * n:].T.reshape(-1, n, 3)
            coll = -np.ones((n, n))
            for i in range(n):
                for j in range(i + 1, n):
                    d = np.linalg.norm(P[:, i] - P[:, j], axis=1)
                    hit = np.where(d < rad[i] + rad[j])[0]
                    if len(hit):
                        coll[i, j] = T[hit[0]]
        else:
            raise ValueError("method must be 'yoshida4' or 'dop853'")
        E = np.array([self.energy(P[k], V[k], pn=pn) for k in range(len(T))])
        scale = abs(E[0]) if E[0] != 0 else 1.0
        collisions = [{"i": i, "j": j, "t": float(coll[i, j]),
                       "names": [self.bodies[i].name, self.bodies[j].name]}
                      for i in range(n) for j in range(i + 1, n) if coll[i, j] >= 0]
        if update:
            self.set_state(P[-1], V[-1])
            self.t += float(T[-1])
        return {"t": T, "pos": P, "vel": V, "energy_drift": (E - E[0]) / scale,
                "collisions": collisions, "method": method, "pn": pn}

    # ---- frames ------------------------------------------------------------
    def transform(self, hist, frame="inertial", pair=None):
        """Re-express a history in another frame.

        'inertial'     as integrated
        'barycentric'  origin at the system's centre of mass
        'corotating'   origin at the barycentre of bodies pair=(i, j), x axis from
                       i to j, z axis along the pair's orbital angular momentum.
                       Velocities are those *seen in the rotating frame*:
                       v_rot = R^T (v - V_origin) - Omega x r_rot.
        """
        P, V = hist["pos"], hist["vel"]
        gm = self.gm
        if frame == "inertial":
            return P.copy(), V.copy()
        if frame == "barycentric":
            r0 = (gm[None, :, None] * P).sum(1) / gm.sum()
            v0 = (gm[None, :, None] * V).sum(1) / gm.sum()
            return P - r0[:, None], V - v0[:, None]
        if frame == "corotating":
            if pair is None:
                raise ValueError("corotating frame needs pair=(i, j)")
            i, j = pair
            mt = gm[i] + gm[j]
            ro = (gm[i] * P[:, i] + gm[j] * P[:, j]) / mt
            vo = (gm[i] * V[:, i] + gm[j] * V[:, j]) / mt
            d, dv = P[:, j] - P[:, i], V[:, j] - V[:, i]
            ex = d / np.linalg.norm(d, axis=1)[:, None]
            h = np.cross(d, dv)
            ez = h / np.linalg.norm(h, axis=1)[:, None]
            ey = np.cross(ez, ex)
            w = np.linalg.norm(h, axis=1) / (d ** 2).sum(1)          # instantaneous spin rate
            R = np.stack([ex, ey, ez], axis=1)                       # rows = new axes
            Pr = np.einsum("fab,fnb->fna", R, P - ro[:, None])
            Vr = np.einsum("fab,fnb->fna", R, V - vo[:, None])
            Vr[:, :, 0] += w[:, None] * Pr[:, :, 1]                  # - Omega x r, Omega = w e_z
            Vr[:, :, 1] -= w[:, None] * Pr[:, :, 0]
            return Pr, Vr
        raise ValueError("unknown frame")

    # ---- units -------------------------------------------------------------
    def unit_system(self, pair=(0, 1)):
        """Non-dimensional units built on a chosen pair, with the conversion
        factors spelled out (this is exactly how CR3BP units are defined)."""
        i, j = pair
        L = float(np.linalg.norm(self.pos[j] - self.pos[i]))
        gm = self.bodies[i].gm + self.bodies[j].gm
        T = float(np.sqrt(L ** 3 / gm))
        return {"length_m": L, "time_s": T, "velocity_m_s": L / T, "gm_m3_s2": gm,
                "mu": min(self.bodies[i].gm, self.bodies[j].gm) / gm,
                "c_nd": self.c / (L / T),
                "explain": [f"1 length unit = separation of {self.bodies[i].name} and {self.bodies[j].name} = {L:.6e} m",
                            f"1 time unit = sqrt(L^3 / G(M1+M2)) = {T:.6e} s  (orbital period / 2pi)",
                            f"1 velocity unit = L/T = {L/T:.6e} m/s",
                            f"speed of light = {self.c/(L/T):.6e} velocity units"]}

    def to_dict(self):
        return {"t": self.t, "bodies": [b.to_dict() for b in self.bodies]}


def two_body(primary="Earth", secondary="Moon", distance=3.844e8):
    """Convenience: a primary with a secondary on a circular orbit, barycentric."""
    s = NBodySystem([PRESETS[primary]])
    s.add_circular_satellite(secondary, 0, distance)
    s.to_barycentric()
    return s
