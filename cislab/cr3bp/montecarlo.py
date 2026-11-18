"""Monte Carlo injection-error experiment (the engine behind the paper).

Question
--------
A spacecraft is placed on the stable manifold of an L2 halo orbit at the point
where that manifold crosses the Moon plane (x = 1 - mu).  With a perfect state
it coasts into the halo for free.  How large may the injection error be before
that stops working?

Experiment
----------
1. Nominal injection states: section crossings of the Moon-bound stable branch
   for several phases around the halo.
2. Add Gaussian errors: position ~ N(0, sigma_r^2) per axis, velocity
   ~ N(0, sigma_v^2) per axis.
3. Coast with NO correction burns and record what happens.
4. Classify under two definitions (both are reported, because the answer
   depends on what you mean by "captured"):

   A  "halo-tracking":  the spacecraft stays inside a tube of radius d_tol
      around the halo for at least n_rev consecutive revolutions.
   B  "Moon-region":    the spacecraft stays within R_B (default: two Hill radii)
      of the Moon, without impact, until the horizon t_entry + n_hold periods
      (t_entry = when the *nominal* trajectory enters the tube; n_hold = n_rev = 1
      by default so that A and B are judged over the same length of time).

   In both, hitting the Moon (or Earth) first is 'impacted'; otherwise 'escaped'.

Common random numbers
---------------------
The SAME unit-Gaussian draws are reused for every error level, every orbit and
every model (Newtonian / 1PN / tighter tolerance).  Differences between two
curves are then caused by the thing we changed, not by sampling luck - a
standard variance-reduction technique.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from .._jit import njit, prange
from . import fastprop as fp
from .fastprop import dop853_step, _rhs, _A, _B, _C, _E3, _E5
from .manifolds import manifold, branch_towards_secondary
from .periodic import PeriodicOrbit
from ..constants import (R_MOON_KM, R_EARTH_KM, LUNIT_EM_KM, C_ND_EM, VUNIT_EM_KMS)

CAPTURED, ESCAPED, IMPACTED = 0, 1, 2
LABELS = {CAPTURED: "captured", ESCAPED: "escaped", IMPACTED: "impacted"}


def hill_radius(mu):
    """Hill radius of the secondary, (mu/3)^(1/3), in units of the primaries' distance."""
    return (mu / 3.0) ** (1.0 / 3.0)


@njit(cache=True)
def _seg_min_dist2(ax, ay, az, bx, by, bz, cx):
    """Squared closest approach to the point (cx,0,0) along the chord a->b."""
    dx, dy, dz = bx - ax, by - ay, bz - az
    L2 = dx * dx + dy * dy + dz * dz
    if L2 == 0.0:
        return (ax - cx) ** 2 + ay * ay + az * az
    s = -((ax - cx) * dx + ay * dy + az * dz) / L2
    s = 0.0 if s < 0.0 else (1.0 if s > 1.0 else s)
    px, py, pz = ax + s * dx - cx, ay + s * dy, az + s * dz
    return px * px + py * py + pz * pz


@njit(cache=True)
def _one(y0, mu, t_max, rtol, atol, hmax, halo, d_tol, R_B, r_moon, r_earth, pn, c,
         A, B, C, E3, E5, out):
    n = 6
    y = y0.copy()
    f = np.empty(n); fnew = np.empty(n); ynew = np.empty(n)
    K = np.empty((13, n))
    _rhs(y, mu, n, pn, c, f)
    h = min(1e-3, hmax)
    t = 0.0
    t_impact, body = -1.0, 0.0
    t_exit = -1.0
    first_entry = -1.0
    longest, run_start, inside = 0.0, 0.0, False
    dmin = 1e60
    d2tol = d_tol * d_tol
    m = halo.shape[0]
    while t < t_max:
        if t + h > t_max:
            h = t_max - t
        err = dop853_step(y, f, h, mu, n, pn, c, rtol, atol, K, ynew, fnew, A, B, C, E3, E5)
        if err > 1.0:
            h *= max(0.2, 0.9 * err ** (-0.125))
            if h < 1e-14:
                t_impact, body = t, 2.0       # step underflow = fell into the Moon's singularity
                break
            continue
        # --- accepted step: test the chord for impacts -----------------------
        if _seg_min_dist2(y[0], y[1], y[2], ynew[0], ynew[1], ynew[2], 1.0 - mu) < r_moon * r_moon:
            t_impact, body = t + h, 2.0
        elif _seg_min_dist2(y[0], y[1], y[2], ynew[0], ynew[1], ynew[2], -mu) < r_earth * r_earth:
            t_impact, body = t + h, 1.0
        t += h
        for i in range(n):
            y[i] = ynew[i]
            f[i] = fnew[i]
        fac = 10.0 if err == 0.0 else min(10.0, 0.9 * err ** (-0.125))
        h = min(h * fac, hmax)
        if t_impact >= 0.0:
            break
        # --- distance to the Moon (definition B) ----------------------------
        r2 = np.sqrt((y[0] - 1.0 + mu) ** 2 + y[1] * y[1] + y[2] * y[2])
        if t_exit < 0.0 and r2 > R_B:
            t_exit = t
        # --- distance to the halo (definition A) ----------------------------
        d2 = 1e300
        for j in range(m):
            q = (y[0] - halo[j, 0]) ** 2 + (y[1] - halo[j, 1]) ** 2 + (y[2] - halo[j, 2]) ** 2
            if q < d2:
                d2 = q
        if d2 < dmin:
            dmin = d2
        if d2 < d2tol:
            if not inside:
                inside, run_start = True, t
                if first_entry < 0.0:
                    first_entry = t
            if t - run_start > longest:
                longest = t - run_start
        else:
            inside = False
        if r2 > 3.0 * R_B and t_exit >= 0.0 and (y[0] * y[0] + y[1] * y[1]) > 6.25:
            break                             # long gone: no need to keep integrating
    out[0], out[1], out[2], out[3], out[4], out[5], out[6] = (t_impact, body, t_exit, first_entry,
                                                              longest, t, np.sqrt(dmin))


@njit(parallel=True, cache=True)
def _batch(Y0, mu, t_max, rtol, atol, hmax, halo, d_tol, R_B, r_moon, r_earth, pn, c,
           A, B, C, E3, E5):
    N = Y0.shape[0]
    out = np.empty((N, 7))
    for i in prange(N):
        _one(Y0[i], mu, t_max, rtol, atol, hmax, halo, d_tol, R_B, r_moon, r_earth, pn, c,
             A, B, C, E3, E5, out[i])
    return out


def fly(Y0, mu, t_max, halo_xyz, d_tol, R_B, rtol=1e-12, atol=1e-14, hmax=0.01, pn=False,
        c=C_ND_EM, r_moon=R_MOON_KM / LUNIT_EM_KM, r_earth=R_EARTH_KM / LUNIT_EM_KM):
    """Propagate many initial states and return the raw outcome table.

    Columns: t_impact (-1 = none), body (1 Earth, 2 Moon), t_exit_B (-1 = never left),
             first_entry into the halo tube (-1 = never), longest dwell in the tube,
             final time, minimum distance to the halo.
    """
    return _batch(np.ascontiguousarray(Y0, dtype=np.float64), float(mu), float(t_max), rtol, atol,
                  hmax, np.ascontiguousarray(halo_xyz, dtype=np.float64), d_tol, R_B, r_moon,
                  r_earth, 1 if pn else 0, float(c), _A, _B, _C, _E3, _E5)


def classify(raw, period, t_horizon, n_rev=1.0, definition="A"):
    """Turn the raw outcome table into labels CAPTURED / ESCAPED / IMPACTED."""
    t_imp, t_exit, dwell = raw[:, 0], raw[:, 2], raw[:, 4]
    lab = np.full(len(raw), ESCAPED)
    if definition == "A":
        ok = dwell >= n_rev * period
        lab[ok] = CAPTURED
        lab[(~ok) & (t_imp >= 0)] = IMPACTED
    elif definition == "B":
        hit = (t_imp >= 0) & (t_imp <= t_horizon)
        left = (t_exit >= 0) & (t_exit <= t_horizon)
        lab[~hit & ~left] = CAPTURED
        lab[hit & ~(left & (t_exit < t_imp))] = IMPACTED
    else:
        raise ValueError("definition must be 'A' or 'B'")
    return lab


@dataclass
class InjectionCase:
    """Everything needed to rerun the experiment for one halo orbit."""
    orbit: PeriodicOrbit
    phases: np.ndarray            # manifold phases used
    states: np.ndarray            # nominal injection states (k, 6) on the Moon plane
    t_flight: np.ndarray          # nominal time from injection to first entry of the tube
    halo_xyz: np.ndarray          # dense samples of the halo (for distance tests)
    d_tol: float
    R_B: float
    t_max: float
    t_horizon: np.ndarray         # per phase: t_flight + n_hold * period
    nominal_dwell: np.ndarray = field(default=None)
    gains: np.ndarray = field(default=None)   # per phase (g_r, g_v): linear error growth
    scan: dict = field(default_factory=dict)  # how many candidate phases were usable


def prepare_case(orbit: PeriodicOrbit, n_phases=8, eps=1e-6, d_tol=5000.0 / LUNIT_EM_KM,
                 R_B=None, n_hold=1.0, n_scan=48, rtol=1e-12, atol=1e-14, hmax=0.01):
    """Find nominal injection states for one orbit and characterise them."""
    mu = orbit.mu
    if R_B is None:
        R_B = 2.0 * hill_radius(mu)
    _, hy = orbit.sample(600)
    halo = hy[:3].T.copy()
    r2max = np.sqrt((halo[:, 0] - 1 + mu) ** 2 + halo[:, 1] ** 2 + halo[:, 2] ** 2).max()
    if r2max + d_tol >= R_B:
        raise ValueError("halo does not fit inside the definition-B region; enlarge R_B")
    sg = branch_towards_secondary(orbit, "stable", eps)
    pts = [p for p in manifold(orbit, "stable", n_scan, eps, sg, section_x=1 - mu, t_max=16.0)
           if p.crossing is not None]
    n_scan_ok = len(pts)
    if not pts:
        raise RuntimeError('no manifold trajectory reaches the Moon plane')
    t_max = float(max(-p.t_cross for p in pts) + 3.0 * orbit.period)
    # Fly every candidate forward once.  A candidate is *usable* only if the
    # unperturbed trajectory really does reach the halo and stays >= 2 periods:
    # trajectories that graze the Moon amplify even round-off error so much that
    # no finite-precision experiment (and no real spacecraft) could use them.
    raw = fly(np.array([p.crossing for p in pts]), mu, t_max, halo, d_tol, R_B, rtol, atol, hmax)
    good = [i for i in range(len(pts)) if raw[i, 3] >= 0 and raw[i, 4] >= 2.0 * orbit.period
            and raw[i, 0] < 0]
    if n_phases is None:                       # use every usable route
        n_phases = len(good)
    if len(good) == 0 or len(good) < n_phases:
        raise RuntimeError(f"only {len(good)} usable manifold trajectories reach the Moon plane")
    pick = [good[k] for k in np.round(np.linspace(0, len(good) - 1, n_phases)).astype(int)]
    pts = [pts[i] for i in pick]
    states = np.array([p.crossing for p in pts])
    raw = raw[pick]
    t_flight = raw[:, 3].copy()
    gains = np.empty((len(pts), 2))
    for i, s in enumerate(states):
        # Linear prediction of error growth from injection to (entry + 1 period):
        # largest singular value of the position rows of the STM, separately for
        # position errors (left 3 columns) and velocity errors (right 3 columns).
        _, Phi = fp.integrate_stm(s, 0.0, t_flight[i] + orbit.period, mu)
        gains[i, 0] = np.linalg.svd(Phi[:3, :3], compute_uv=False)[0]
        gains[i, 1] = np.linalg.svd(Phi[:3, 3:], compute_uv=False)[0]
    case = InjectionCase(orbit, np.array([p.phase for p in pts]), states, t_flight, halo, d_tol, R_B,
                         t_max, t_flight + n_hold * orbit.period, raw[:, 4].copy(), gains)
    case.scan = {"n_scan": n_scan, "n_reach_section": n_scan_ok, "n_usable": len(good)}
    return case


def linear_sigma50(gains, d_tol, sigmas=None):
    """Pooled 50% tolerance predicted by linear theory for a set of routes.

    For one route with growth factor g, a Gaussian error of size sigma is captured
    with probability erf(d_tol / (g sigma sqrt 2)).  Averaging that over routes
    gives the pooled capture curve; we return the sigma where it crosses 1/2.
    """
    from scipy.special import erf
    g = np.asarray(gains, float)
    if sigmas is None:
        sigmas = np.logspace(-11, -1, 400)
    p = np.array([np.mean(erf(d_tol / (g * s * np.sqrt(2)))) for s in sigmas])
    return sigma50(sigmas, p)


def unit_draws(n_samples, seed=20261010):
    """The common random numbers: n_samples x 6 standard normal values."""
    return np.random.default_rng(seed).standard_normal((n_samples, 6))


def run_levels(case: InjectionCase, sigma_r, sigma_v, Z, pn=False, rtol=1e-12, atol=1e-14,
               hmax=0.01, n_rev=1.0):
    """Capture statistics for paired lists of (sigma_r, sigma_v).

    Returns dict with arrays shaped (levels, phases, samples) of labels under
    definitions A and B.
    """
    sigma_r, sigma_v = np.atleast_1d(sigma_r).astype(float), np.atleast_1d(sigma_v).astype(float)
    L, P, S = len(sigma_r), len(case.states), len(Z)
    Y0 = np.empty((L, P, S, 6))
    for l in range(L):
        Y0[l, :, :, :3] = case.states[:, None, :3] + sigma_r[l] * Z[None, :, :3]
        Y0[l, :, :, 3:] = case.states[:, None, 3:] + sigma_v[l] * Z[None, :, 3:]
    raw = fly(Y0.reshape(-1, 6), case.orbit.mu, case.t_max, case.halo_xyz, case.d_tol, case.R_B,
              rtol, atol, hmax, pn).reshape(L, P, S, 7)
    labA = np.empty((L, P, S), dtype=np.int8)
    labB = np.empty((L, P, S), dtype=np.int8)
    for p in range(P):
        r = raw[:, p].reshape(-1, 7)
        labA[:, p] = classify(r, case.orbit.period, case.t_horizon[p], n_rev, "A").reshape(L, S)
        labB[:, p] = classify(r, case.orbit.period, case.t_horizon[p], n_rev, "B").reshape(L, S)
    return {"A": labA, "B": labB, "raw": raw.astype(np.float32)}


def capture_probability(labels):
    """Fraction captured at each level, pooled over phases and samples."""
    return (labels == CAPTURED).mean(axis=(1, 2))


def sigma50(sigmas, prob, level=0.5):
    """Error size at which capture probability falls through ``level``
    (log-linear interpolation).  NaN if the curve never crosses."""
    s, p = np.asarray(sigmas, float), np.asarray(prob, float)
    for i in range(len(s) - 1):
        if p[i] >= level > p[i + 1]:
            w = (p[i] - level) / (p[i] - p[i + 1])
            return float(np.exp(np.log(s[i]) + w * (np.log(s[i + 1]) - np.log(s[i]))))
    return float("nan")


def sigma50_bootstrap(sigmas, labels, n_boot=400, seed=7, level=0.5):
    """Bootstrap 16-84% interval of sigma50, resampling the random draws (kept
    common across levels, exactly as in the experiment) and the phases."""
    rng = np.random.default_rng(seed)
    cap = labels == CAPTURED
    L, P, S = cap.shape
    vals = []
    for _ in range(n_boot):
        ip = rng.integers(0, P, P)
        is_ = rng.integers(0, S, S)
        vals.append(sigma50(sigmas, cap[:, ip][:, :, is_].mean(axis=(1, 2)), level))
    vals = np.array(vals)
    vals = vals[np.isfinite(vals)]
    if len(vals) < n_boot // 2:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 16)), float(np.percentile(vals, 84))


def to_km(x_nd):
    return x_nd * LUNIT_EM_KM


def to_m_per_s(v_nd):
    return v_nd * VUNIT_EM_KMS * 1000.0
