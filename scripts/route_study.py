"""Route study: how much do pooled tolerances depend on WHICH manifold routes are sampled?

    python scripts/route_study.py        # ~20 minutes on 2 cores

Motivation.  The robustness run with a ten-times smaller manifold step (eps)
changed the pooled 50% tolerance by more than its bootstrap range.  A smaller
eps does not move the manifold tube, but it does pick 12 *different*
trajectories on it, and tolerance varies by orders of magnitude from route to
route.  So 12 routes are too few to pin down a pooled value.  This script
  (1) uses the linear (STM) prediction, validated route-by-route by the Monte
      Carlo, on EVERY usable route of a 96-trajectory scan, for 24 halos; and
  (2) repeats the Monte Carlo itself with every usable route of a 48-scan, for
      three halos and for both eps values, to test (1) and the eps-independence.
Writes data/route_study.json.
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cislab.constants import MU_EARTH_MOON as MU, LUNIT_EM_KM as KM, VUNIT_EM_KMS
from cislab.cr3bp.families import build_family, member_with_z
from cislab.cr3bp import montecarlo as mc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MMS = VUNIT_EM_KMS * 1e6
t0 = time.time()
fam = build_family(MU)
out = {"linear": [], "dense_mc": []}
for z in np.round(np.linspace(0.012, 0.196, 24), 4):
    orb = member_with_z(fam, float(z))
    row = {"z_far": float(z), "Az_km": orb.amplitudes()[2] * KM, "stability": orb.stability_index}
    for eps in (1e-6, 1e-7):
        try:
            c = mc.prepare_case(orb, n_phases=None, n_scan=96, eps=eps)
        except Exception as e:
            row[f"error_{eps:g}"] = str(e)
            continue
        tag = "" if eps == 1e-6 else "_eps7"
        tol_r = c.d_tol / (0.6745 * c.gains[:, 0]) * KM * 1e3            # per-route, metres
        row.update({f"n_usable{tag}": int(c.scan["n_usable"]), f"n_scan{tag}": int(c.scan["n_scan"]),
                    f"n_reach{tag}": int(c.scan["n_reach_section"]),
                    f"pos_m{tag}": mc.linear_sigma50(c.gains[:, 0], c.d_tol) * KM * 1e3,
                    f"vel_mms{tag}": mc.linear_sigma50(c.gains[:, 1], c.d_tol) * MMS,
                    f"route_q10_m{tag}": float(np.percentile(tol_r, 10)), f"route_q50_m{tag}": float(np.percentile(tol_r, 50)),
                    f"route_q90_m{tag}": float(np.percentile(tol_r, 90)), f"route_best_m{tag}": float(tol_r.max()),
                    f"route_worst_m{tag}": float(tol_r.min())})
    out["linear"].append(row)
    print(f"z={z:.3f} k={row['stability']:7.1f} linear pooled pos: {row.get('pos_m', float('nan')):9.1f} m (eps/10: {row.get('pos_m_eps7', float('nan')):9.1f})  "
          f"usable {row.get('n_usable')}/{row.get('n_scan')}", flush=True)

SIG = np.logspace(-10, -3, 22)
Z = mc.unit_draws(80, 20261010)
for z in (0.04, 0.10, 0.16):
    orb = member_with_z(fam, z)
    for eps in (1e-6, 1e-7):
        c = mc.prepare_case(orb, n_phases=None, n_scan=48, eps=eps)
        r = mc.run_levels(c, SIG, 0 * SIG, Z)
        s = mc.sigma50(SIG, mc.capture_probability(r["A"]))
        lo, hi = mc.sigma50_bootstrap(SIG, r["A"])
        out["dense_mc"].append({"z_far": z, "eps": eps, "n_routes": len(c.states), "n_traj": int(r["A"].size),
                                "pos_m": s * KM * 1e3, "lo_m": lo * KM * 1e3, "hi_m": hi * KM * 1e3,
                                "linear_m": mc.linear_sigma50(c.gains[:, 0], c.d_tol) * KM * 1e3})
        print("dense MC", out["dense_mc"][-1], flush=True)
out["seconds"] = time.time() - t0
json.dump(out, open(os.path.join(ROOT, "data", "route_study.json"), "w"), indent=1)
print(f"done in {time.time()-t0:.0f}s")
