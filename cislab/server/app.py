"""Local server: the thin layer between the physics core and the browser.

    python -m cislab.server            # then open http://127.0.0.1:8000

Rules this file follows:
* It contains NO physics.  Every number comes from a function in ``cislab``.
* Every response that involves integration carries an accuracy readout
  (Jacobi drift for the CR3BP, energy drift for the N-body problem).
"""
from __future__ import annotations
import glob
import json
import os
from functools import lru_cache
from typing import Optional

import numpy as np
from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .. import constants as K
from ..constants import Body, PRESETS
from ..cr3bp import (lagrange_points, classify_point, linear_frequencies, jacobi, omega)
from ..cr3bp import fastprop as fp
from ..cr3bp.hill import critical_jacobi, hill_slice, regime as hill_regime
from ..cr3bp.families import build_family, member_with_z, family_table
from ..cr3bp.richardson import Richardson
from ..cr3bp.manifolds import (manifold, branch_towards_secondary, section_curve, convergence_check,
                               eigen_directions)
from ..cr3bp.poincare import planar_map
from ..cr3bp import montecarlo as mc
from ..nbody.system import NBodySystem
from ..relativity import regime as rg, clocks, light, schwarzschild as sw
from ..relativity.eih import newton_vs_1pn
from .. import fields

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WEB = os.path.join(ROOT, "web")
app = FastAPI(title="Cislunar Lab", version="1.0")


def _clean(o):
    """Make NumPy objects JSON-serialisable (NaN/inf become null)."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if np.isfinite(f) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def _r(a, nd=7):
    """Round coordinates sent for drawing (keeps responses small)."""
    return np.round(np.asarray(a, float), nd)


# --------------------------------------------------------------------------- meta
@app.get("/api/meta")
def meta():
    return _clean({
        "mu_earth_moon": K.MU_EARTH_MOON, "lunit_km": K.LUNIT_EM_KM, "tunit_s": K.TUNIT_EM_S,
        "vunit_kms": K.VUNIT_EM_KMS, "c": K.C, "G": K.G, "c_nd": K.C_ND_EM,
        "r_moon": K.R_MOON_KM / K.LUNIT_EM_KM, "r_earth": K.R_EARTH_KM / K.LUNIT_EM_KM,
        "presets": {k: v.to_dict() for k, v in PRESETS.items()}, "sources": K.SOURCES,
    })


# --------------------------------------------------------------------------- CR3BP
@app.get("/api/cr3bp/lagrange")
def api_lagrange(mu: float = K.MU_EARTH_MOON):
    if not 0 < mu <= 0.5:
        raise HTTPException(400, "mu must be in (0, 0.5]")
    out = {}
    for name in ("L1", "L2", "L3", "L4", "L5"):
        d = classify_point(mu, name)
        if name in ("L1", "L2", "L3"):
            lam, nu, wz = linear_frequencies(mu, name)
            d.update({"lambda": lam, "nu": nu, "omega_z": wz,
                      "e_fold_days": K.TUNIT_EM_S / lam / 86400 if abs(mu - K.MU_EARTH_MOON) < 1e-6 else None})
        out[name] = d
    return _clean({"mu": mu, "points": out, "critical_jacobi": critical_jacobi(mu)})


@app.get("/api/cr3bp/hill")
def api_hill(mu: float = K.MU_EARTH_MOON, C: float = 3.15, n: int = 220, lim: float = 1.5,
             zlevels: int = 0):
    """Zero-velocity curve data: the plane grid for tinting, plus contour polylines
    (also at several heights z, which together sketch the 3D zero-velocity surface)."""
    import contourpy
    u, v, W = hill_slice(mu, C, "xy", 0.0, (-lim, lim, -lim, lim), n)
    zs = [0.0] + [float(z) for z in np.linspace(0.08, 0.9, zlevels)] if zlevels else [0.0]
    lines = []
    for z in zs:
        _, _, Wz = hill_slice(mu, C, "xy", z, (-lim, lim, -lim, lim), 200)
        uu = np.linspace(-lim, lim, 200)
        for seg in contourpy.contour_generator(uu, uu, Wz).lines(0.0):
            if len(seg) > 3:
                lines.append({"z": z, "pts": _r(seg, 4)})
    # x-z slice too: shows how the surface closes over the poles
    uu = np.linspace(-lim, lim, 200)
    _, _, Wxz = hill_slice(mu, C, "xz", 0.0, (-lim, lim, -lim, lim), 200)
    xz = [{"pts": _r(seg, 4)} for seg in contourpy.contour_generator(uu, uu, Wxz).lines(0.0) if len(seg) > 3]
    return _clean({"mu": mu, "C": C, "lim": lim, "n": n, "grid": _r(np.clip(W, -1, 1), 3).ravel(),
                   "contours": lines, "contours_xz": xz, "regime": hill_regime(mu, C),
                   "critical": critical_jacobi(mu)})


@lru_cache(maxsize=4)
def _family(mu: float, point: str):
    return build_family(mu, point)


def _orbit_payload(o, n=240):
    t, y = o.sample(n)
    d = o.to_dict()
    d.update({"points": _r(y.T, 6), "multipliers": [[float(m.real), float(m.imag)] for m in o.multipliers],
              "Az_km": d["Az"] * K.LUNIT_EM_KM, "period_days": o.period * K.TUNIT_EM_S / 86400,
              "jacobi_jpl": o.jacobi - o.mu * (1 - o.mu)})
    return d


@app.get("/api/cr3bp/family")
def api_family(mu: float = K.MU_EARTH_MOON, point: str = "L2", every: int = 6):
    fam = _family(mu, point)
    tab = family_table(fam)
    ghosts = [_r(o.sample(90)[1][:3].T, 5) for o in fam[::every]]
    z = np.abs(tab[:, 1])
    return _clean({"mu": mu, "point": point, "table": tab, "columns": ["x0", "z0", "vy0", "period", "jacobi", "stability"],
                   "ghosts": ghosts, "z_max": float(z.max()), "z_min": float(z.min()),
                   "i_peak": int(np.argmax(z))})


@app.get("/api/cr3bp/halo")
def api_halo(z: float, mu: float = K.MU_EARTH_MOON, point: str = "L2", branch: str = "first"):
    try:
        o = member_with_z(_family(mu, point), z, branch)
    except ValueError as e:
        raise HTTPException(400, str(e))
    d = _orbit_payload(o)
    # Richardson's analytic guess for the same out-of-plane amplitude, for comparison
    try:
        R = Richardson(mu, point)
        az = d["Az"]
        ts = np.linspace(0, R.period(az / R.gamma), 160)
        d["richardson"] = _r(np.array([R.state(az, t, north=(point != "L2"))[:3] for t in ts]), 6)
    except Exception:
        d["richardson"] = None
    return _clean(d)


@app.get("/api/cr3bp/manifold")
def api_manifold(z: float, kind: str = "stable", towards_moon: bool = True, eps: float = 1e-6,
                 n: int = 40, mu: float = K.MU_EARTH_MOON, t_max: float = 8.0, cut: bool = True,
                 branch: str = "first"):
    o = member_with_z(_family(mu, "L2"), z, branch)
    try:
        eigen_directions(o)
    except ValueError as e:
        raise HTTPException(400, str(e))
    sg = branch_towards_secondary(o, kind, eps)
    if not towards_moon:
        sg = -sg
    pts = manifold(o, kind, n, eps, sg, t_max=t_max, section_x=(1 - mu) if (cut and towards_moon) else None)
    tubes = []
    for p in pts:
        step = max(1, len(p.t) // 160)
        tubes.append({"phase": p.phase, "impact": p.impact, "pts": _r(p.y[::step, :3], 5)})
    ph, cur = section_curve(pts)
    conv = convergence_check(o, kind, eps, sg) if (cut and towards_moon) else None
    return _clean({"kind": kind, "eps": eps, "eps_km": eps * K.LUNIT_EM_KM, "tubes": tubes,
                   "section": {"phase": ph, "states": _r(cur, 6)} if len(cur) else None,
                   "convergence": conv, "stability": o.stability_index})


class PropagateReq(BaseModel):
    state: list[float]
    t: float
    mu: float = K.MU_EARTH_MOON
    pn: bool = False
    compare_pn: bool = False


@app.post("/api/cr3bp/propagate")
def api_propagate(req: PropagateReq):
    s0 = np.array(req.state, float)
    T, Y = fp.trajectory(s0, 0.0, req.t, req.mu, hmax=0.01, pn=req.pn, c=K.C_ND_EM)
    Cj = jacobi(Y.T, req.mu)
    out = {"t": _r(T, 6), "states": _r(Y, 12), "jacobi0": float(Cj[0]),
           "jacobi_drift": float(np.max(np.abs(Cj - Cj[0])) / abs(Cj[0])),
           "note": "With the 1PN term on, the Newtonian Jacobi constant is no longer exactly "
                   "conserved, so its 'drift' then measures physics, not numerical error." if req.pn else ""}
    if req.compare_pn:
        y = s0.copy()
        gap = [0.0]
        for k in range(1, len(T)):
            y = fp.integrate(y, T[k - 1], T[k], req.mu, pn=True, c=K.C_ND_EM)
            gap.append(float(np.linalg.norm(y[:3] - Y[k, :3])))
        out["pn_gap_m"] = _r(np.array(gap) * K.LUNIT_EM_KM * 1000, 6)
        ax = np.array(fp.pn_test_particle(*s0, req.mu, K.C_ND_EM))
        from ..cr3bp.core import grad_omega
        out["pn_ratio"] = float(np.linalg.norm(ax) / np.linalg.norm(grad_omega(*s0[:3], req.mu) - [s0[0], s0[1], 0]))
    return _clean(out)


@app.get("/api/cr3bp/poincare")
def api_poincare(C: float = 3.15, mu: float = K.MU_EARTH_MOON, n: int = 28, t_max: float = 40.0,
                 xmin: float = 0.86, xmax: float = 1.14):
    xs = np.linspace(xmin, xmax, n)
    xs = xs[np.abs(xs - (1 - mu)) > 0.01]
    orbits = planar_map(mu, C, xs, t_max=t_max)
    return _clean({"C": C, "orbits": [{"x0": x0, "pts": _r(p, 5)} for x0, p in orbits]})


class MCReq(BaseModel):
    z: float = 0.10
    sigma_r_km: float = 1.0
    sigma_v_mms: float = 0.0
    n: int = 120
    pn: bool = False
    n_draw: int = 36


@lru_cache(maxsize=8)
def _case(z: float):
    return mc.prepare_case(member_with_z(_family(K.MU_EARTH_MOON, "L2"), z), n_phases=6)


@app.post("/api/cr3bp/montecarlo")
def api_mc(req: MCReq):
    """A small live Monte Carlo for the interface (the paper uses the big batch run)."""
    case = _case(round(req.z, 4))
    sr = req.sigma_r_km / K.LUNIT_EM_KM
    sv = req.sigma_v_mms / 1e6 / K.VUNIT_EM_KMS
    Z = mc.unit_draws(max(10, min(req.n, 600)))
    res = mc.run_levels(case, [sr], [sv], Z, pn=req.pn)
    A, B = res["A"][0], res["B"][0]
    draw = []
    k = 0
    P, S = A.shape
    for i in range(min(req.n_draw, P * S)):
        p, s = i % P, (i * 7) % S
        y0 = case.states[p].copy()
        y0[:3] += sr * Z[s, :3]
        y0[3:] += sv * Z[s, 3:]
        T, Y = fp.trajectory(y0, 0.0, float(case.t_horizon[p] + 0.5 * case.orbit.period), K.MU_EARTH_MOON, hmax=0.03,
                             rtol=1e-10, atol=1e-12)
        r2 = np.sqrt((Y[:, 0] - 1 + K.MU_EARTH_MOON) ** 2 + Y[:, 1] ** 2 + Y[:, 2] ** 2)
        hit = np.where(r2 < K.R_MOON_KM / K.LUNIT_EM_KM)[0]
        far = np.where(np.abs(Y[:, 0]) + np.abs(Y[:, 1]) > 3)[0]
        end = min([len(Y)] + ([hit[0] + 1] if len(hit) else []) + ([far[0]] if len(far) else []))
        draw.append({"A": int(A[p, s]), "B": int(B[p, s]), "pts": _r(Y[:end:2, :3], 5)})
    frac = lambda L: {mc.LABELS[k]: float((L == k).mean()) for k in (0, 1, 2)}
    return _clean({"A": frac(A), "B": frac(B), "n": int(A.size), "tracks": draw,
                   "d_tol_km": case.d_tol * K.LUNIT_EM_KM, "R_B_km": case.R_B * K.LUNIT_EM_KM,
                   "injection": _r(case.states[:, :3], 6), "stability": case.orbit.stability_index,
                   "standard_error": float(np.sqrt(0.25 / A.size))})


@app.get("/api/cr3bp/paper_results")
def api_paper_results():
    """Summary of the batch Monte Carlo (written by scripts/analyze_experiment.py)."""
    p = os.path.join(ROOT, "data", "mc_summary.json")
    if not os.path.exists(p):
        raise HTTPException(404, "Run scripts/run_experiment.py and scripts/analyze_experiment.py first")
    return json.load(open(p))


# --------------------------------------------------------------------------- N-body
class BodyIn(BaseModel):
    name: str
    gm: float
    radius: float
    pos: list[float]
    vel: list[float]
    spin: list[float] = [0, 0, 0]


class NBodyReq(BaseModel):
    bodies: list[BodyIn]
    t_end: float = 0.0
    method: str = "dop853"
    dt: Optional[float] = None
    pn: bool = False
    compare: bool = False
    frame: str = "inertial"
    pair: Optional[list[int]] = None
    n_frames: int = 400


def _system(bodies):
    return NBodySystem([Body(b.name, b.gm, b.radius, tuple(b.pos), tuple(b.vel), tuple(b.spin)) for b in bodies])


def _analysis(s: NBodySystem, pair=None):
    names = [b.name for b in s.bodies]
    r0, v0 = s.barycenter()
    rep = rg.regime_report(s.gm, s.pos - r0, s.vel - v0, s.radius, names)
    rates = clocks.body_clock_rates(s.gm, s.pos, s.vel - v0, s.radius)
    for r, cr, b in zip(rep, rates, s.bodies):
        r["clock_rate_offset"] = float(cr)
        r["clock_us_per_day"] = float(cr * 86400e6)
        r["rs"] = b.schwarzschild_radius
        w = float(np.linalg.norm(b.spin))
        r["spin_period_s"] = 2 * np.pi / w if w > 0 else None
        # Rotation numbers are *estimates* shown for context only (see docs):
        r["surface_speed_over_c"] = w * b.radius / s.c
        r["breakup_fraction"] = w / np.sqrt(b.gm / b.radius ** 3) if b.radius > 0 else None
    out = {"regime": rep, "energy_over_G": s.energy(), "momentum": s.momentum(),
           "angular_momentum": s.angular_momentum()}
    if pair and len(s.bodies) >= 2:
        out["pair"] = fields.pair_report(s.bodies, tuple(pair))
        out["units"] = s.unit_system(tuple(pair))
    return out


@app.post("/api/nbody/analyze")
def api_analyze(req: NBodyReq):
    return _clean(_analysis(_system(req.bodies), req.pair))


@app.post("/api/nbody/run")
def api_nbody_run(req: NBodyReq):
    s = _system(req.bodies)
    if len(s.bodies) < 1:
        raise HTTPException(400, "add at least one body")
    try:
        if req.compare:
            cmp_ = newton_vs_1pn(s, req.t_end, n_frames=req.n_frames)
            h = cmp_["pn"] if req.pn else cmp_["newton"]
            other = cmp_["newton"] if req.pn else cmp_["pn"]
        else:
            h = s.run(req.t_end, dt=req.dt, method=req.method, pn=req.pn, n_frames=req.n_frames, update=False)
            cmp_, other = None, None
    except (ValueError, RuntimeError) as e:
        raise HTTPException(400, str(e))
    pair = tuple(req.pair) if req.pair else None
    P, V = s.transform(h, req.frame, pair)
    out = {"t": h["t"], "pos": P, "vel": V, "energy_drift": h["energy_drift"],
           "max_energy_drift": float(np.max(np.abs(h["energy_drift"]))), "collisions": h["collisions"],
           "method": h["method"], "pn": h["pn"], "frame": req.frame,
           "final": {"pos": h["pos"][-1], "vel": h["vel"][-1]}}
    if cmp_ is not None:
        Po, _ = s.transform(other, req.frame, pair)
        out["other_pos"] = Po
        out["divergence_m"] = cmp_["divergence"]
    return _clean(out)


class FieldReq(BaseModel):
    bodies: list[BodyIn]
    center: list[float] = [0, 0, 0]
    half_width: float = 1.0
    n: int = 72
    probe: Optional[list[float]] = None
    n_lines: int = 0


@app.post("/api/nbody/field")
def api_field(req: FieldReq):
    s = _system(req.bodies)
    sheet = fields.potential_sheet(s.gm, s.pos, req.center, req.half_width, req.n)
    out = {"u": sheet["u"], "phi": sheet["phi"].ravel(), "floor": sheet["floor"], "n": req.n,
           "phi_max": float(sheet["phi"].max())}
    if req.n_lines:
        ang = np.linspace(0, 2 * np.pi, req.n_lines, endpoint=False)
        c = np.array(req.center)
        seeds = [c + 0.97 * req.half_width * np.array([np.cos(a), np.sin(a), 0]) for a in ang]
        out["lines"] = [_r(l[::3], 9 - int(np.log10(req.half_width))) for l in
                        fields.field_lines(s.gm, s.pos, seeds, req.half_width / 150, 600, s.radius)]
    if req.probe is not None:
        w, V = fields.tidal_principal_axes(req.probe, s.gm, s.pos)
        out["tidal"] = {"values": w, "axes": V.T, "g": fields.field(np.array([req.probe]), s.gm, s.pos)[0],
                        "phi": float(fields.potential(np.array([req.probe]), s.gm, s.pos)[0]),
                        "eps": float(clocks.potential_over_c2(req.probe, s.gm, s.pos))}
    return _clean(out)


class RayReq(BaseModel):
    bodies: list[BodyIn]
    origin: list[float]
    direction: list[float]
    spread: float = 0.0
    n_rays: int = 9
    length: float = 1.0


@app.post("/api/nbody/rays")
def api_rays(req: RayReq):
    s = _system(req.bodies)
    d = np.array(req.direction, float)
    d /= np.linalg.norm(d)
    side = np.cross([0, 0, 1.0], d)
    if np.linalg.norm(side) < 1e-9:
        side = np.array([1.0, 0, 0])
    side /= np.linalg.norm(side)
    rays = []
    for k in np.linspace(-1, 1, req.n_rays):
        o = np.array(req.origin, float) + k * req.spread * side
        r = light.trace_ray(o, d, s.gm, s.pos, req.length, n=500, radius=s.radius)
        turn = float(np.arccos(np.clip(r["direction"] @ d, -1, 1)))
        rays.append({"pts": r["points"][::4], "max_eps": r["max_eps"], "blocked": r["blocked"],
                     "deflection_arcsec": turn / K.ARCSEC, "trustworthy": bool(r["max_eps"] < 1e-2)})
    return _clean({"rays": rays})


# --------------------------------------------------------------------------- Schwarzschild
@app.get("/api/schwarzschild/orbit")
def api_sch_orbit(rp: float = 8.0, ra: float = 20.0, n_orbits: float = 4.0):
    try:
        g = sw.bound_orbit(rp, ra, n_orbits=n_orbits, n=1800)
    except ValueError as e:
        raise HTTPException(400, str(e))
    a = 0.5 * (rp + ra)
    e = (ra - rp) / (ra + rp)
    return _clean({"x": _r(g["x"], 4), "y": _r(g["y"], 4), "tau": _r(g["tau"], 3), "t": _r(g["t"], 3),
                   "r": _r(g["r"], 4), "plunged": g["plunged"], "E": g["E"], "L": g["L"],
                   "precession_deg_per_orbit": np.degrees(g["precession_per_orbit"]),
                   "weak_field_deg_per_orbit": np.degrees(6 * np.pi / (a * (1 - e * e))),
                   "isco": sw.isco_numeric(), "photon_sphere": sw.photon_sphere_numeric(), "horizon": 2.0})


@app.get("/api/schwarzschild/rays")
def api_sch_rays(bmin: float = 3.0, bmax: float = 14.0, n: int = 12, x0: float = -40.0):
    rays = []
    for b in np.linspace(bmin, bmax, n):
        r = sw.light_ray(x0, float(b), 0.0)
        d = sw.light_deflection(float(b))
        rays.append({"b": float(b), "x": _r(r["x"][::3], 3), "y": _r(r["y"][::3], 3), "captured": r["captured"],
                     "deflection_deg": np.degrees(d) if np.isfinite(d) else None,
                     "weak_deg": np.degrees(sw.light_deflection_weak(float(b)))})
    return _clean({"rays": rays, "b_crit": 3 * np.sqrt(3)})


@app.get("/api/checks")
def api_checks():
    """Headline verification numbers, computed live so the interface can show them."""
    from ..relativity.eih import mercury_precession_arcsec_per_century
    M = K.GM_SUN / K.C ** 2
    return _clean({
        "gps_us_per_day": clocks.gps_clock_offset()["microseconds_per_day"],
        "mercury_arcsec_per_century": mercury_precession_arcsec_per_century(),
        "sun_limb_deflection_arcsec": sw.light_deflection(6.957e8 / M) / K.ARCSEC,
        "isco_over_M": sw.isco_numeric(),
    })


# --------------------------------------------------------------------------- live N-body stream
@app.websocket("/ws/nbody")
async def ws_nbody(ws: WebSocket):
    """Live session.  The browser keeps a buffer of frames and asks for the next
    chunk when playback nears its end; edits simply replace the system state.

    messages in:  {"cmd": "set", "bodies": [...]}
                  {"cmd": "advance", "t": seconds (may be negative), "n_frames": 60,
                   "method": "dop853"|"yoshida4", "dt": .., "pn": bool}
    messages out: {"ok": true, "t0":.., "t": [...], "pos": [...], "vel": [...],
                   "energy_drift": .., "collisions": [...]}
    """
    await ws.accept()
    system = NBodySystem()
    e_ref = None
    try:
        while True:
            msg = json.loads(await ws.receive_text())
            try:
                if msg["cmd"] == "set":
                    system = _system([BodyIn(**b) for b in msg["bodies"]])
                    system.t = float(msg.get("t", 0.0))
                    e_ref = None
                    await ws.send_text(json.dumps({"ok": True, "set": len(system.bodies)}))
                elif msg["cmd"] == "advance":
                    pn = bool(msg.get("pn", False))
                    method = "dop853" if pn else msg.get("method", "dop853")
                    t0 = system.t
                    if e_ref is None:
                        e_ref = system.energy(pn=pn)
                    h = system.run(float(msg["t"]), dt=msg.get("dt"), method=method, pn=pn,
                                   n_frames=int(msg.get("n_frames", 60)))
                    e_now = system.energy(pn=pn)
                    await ws.send_text(json.dumps(_clean({
                        "ok": True, "t0": t0, "t": h["t"] + t0, "pos": h["pos"], "vel": h["vel"],
                        "energy_drift": (e_now - e_ref) / (abs(e_ref) or 1.0),
                        "collisions": h["collisions"], "method": method, "pn": pn})))
                else:
                    await ws.send_text(json.dumps({"ok": False, "error": "unknown cmd"}))
            except Exception as e:                       # report, keep the session alive
                await ws.send_text(json.dumps({"ok": False, "error": str(e)}))
    except WebSocketDisconnect:
        return


# --------------------------------------------------------------------------- static front end
@app.get("/")
def index():
    return FileResponse(os.path.join(WEB, "index.html"))


app.mount("/", StaticFiles(directory=WEB), name="web")
