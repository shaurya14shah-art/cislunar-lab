"""Run the full Monte Carlo injection-error experiment and save raw results.

    python scripts/run_experiment.py            # full run (about an hour on 2 cores)
    python scripts/run_experiment.py --quick    # small run for a smoke test

Results are written one orbit at a time to data/mc/*.npz, so an interrupted run
resumes where it stopped.  Seeds are fixed: rerunning reproduces every number.
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cislab.constants import MU_EARTH_MOON as MU, LUNIT_EM_KM
from cislab.cr3bp.families import build_family, member_with_z
from cislab.cr3bp import montecarlo as mc
from cislab.cr3bp import fastprop as fp
from cislab.constants import C_ND_EM

ap = argparse.ArgumentParser()
ap.add_argument("--quick", action="store_true")
ap.add_argument("--out", default="data/mc")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

Z_FAR = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.195]
# Tolerance varies far more between routes than between random draws, so the budget
# goes to routes: 24 routes x 150 draws per error level.
N_PHASES, N_SAMPLES = (4, 60) if args.quick else (24, 150)
SIG = np.logspace(-10, -3, 8 if args.quick else 22)       # non-dimensional error levels
SEED = 20261010
Z = mc.unit_draws(N_SAMPLES, SEED)
if args.quick:
    Z_FAR = [0.06, 0.14]

t0 = time.time()
family = build_family(MU)
print(f"family of {len(family)} orbits built in {time.time()-t0:.1f}s", flush=True)


def pn_nominal_shift(case):
    """Deterministic effect of 1PN on the nominal coast: largest position gap
    between the Newtonian and the 1PN trajectory from the same injection state,
    up to (entry + 1 period)."""
    out = []
    for s, tf in zip(case.states, case.t_flight):
        tend = tf + case.orbit.period
        Tn, Yn = fp.trajectory(s, 0, tend, MU, hmax=0.01)
        gap = 0.0
        y = s.copy()
        for k in range(1, len(Tn)):
            y = fp.integrate(y, Tn[k - 1], Tn[k], MU, pn=True, c=C_ND_EM)
            gap = max(gap, float(np.linalg.norm(y[:3] - Yn[k, :3])))
        out.append(gap)
    return np.array(out)


def save(name, case, **arrs):
    np.savez_compressed(os.path.join(args.out, name), sig=SIG, phases=case.phases,
                        states=case.states, t_flight=case.t_flight, t_horizon=case.t_horizon,
                        gains=case.gains, period=case.orbit.period, state0=case.orbit.state0,
                        stability=case.orbit.stability_index, jacobi=case.orbit.jacobi,
                        d_tol=case.d_tol, R_B=case.R_B, nominal_dwell=case.nominal_dwell,
                        scan=np.array([case.scan['n_scan'], case.scan['n_reach_section'], case.scan['n_usable']]), **arrs)


for z in Z_FAR:
    name = f"halo_z{z:.3f}.npz"
    if os.path.exists(os.path.join(args.out, name)):
        print("skip", name, flush=True)
        continue
    t1 = time.time()
    try:
        orb = member_with_z(family, z)
        case = mc.prepare_case(orb, n_phases=N_PHASES)
    except Exception as e:                 # report and move on; never hide a failure
        print(f"z_far={z:.3f} SKIPPED: {e}", flush=True)
        continue
    zero = np.zeros_like(SIG)
    out = {}
    for tag, sr, sv, pn in (("pos", SIG, zero, False), ("vel", zero, SIG, False),
                            ("pos_pn", SIG, zero, True), ("vel_pn", zero, SIG, True)):
        r = mc.run_levels(case, sr, sv, Z, pn=pn)
        out[tag + "_raw"] = r["raw"]
    out["pn_shift"] = pn_nominal_shift(case)
    amp = orb.amplitudes()
    save(name, case, amp=np.array(amp), **out)
    print(f"z_far={z:.3f}  k={orb.stability_index:8.2f}  T={orb.period:.4f}  "
          f"done in {time.time()-t1:.0f}s", flush=True)

# ---- robustness runs on three representative orbits -------------------------
for z in ([0.06] if args.quick else [0.04, 0.10, 0.16]):
    name = f"robust_z{z:.3f}.npz"
    if os.path.exists(os.path.join(args.out, name)):
        continue
    t1 = time.time()
    orb = member_with_z(family, z)
    zero = np.zeros_like(SIG)
    out = {}
    case = mc.prepare_case(orb, n_phases=N_PHASES)
    # (a) every integration tolerance halved, step cap halved
    out["tight_pos_raw"] = mc.run_levels(case, SIG, zero, Z, rtol=5e-13, atol=5e-15, hmax=0.005)["raw"]
    # (b) capture tube half / double as wide
    for fac in (0.5, 2.0):
        c2 = mc.prepare_case(orb, n_phases=N_PHASES, d_tol=fac * case.d_tol)
        out[f"dtol{fac}_pos_raw"] = mc.run_levels(c2, SIG, zero, Z)["raw"]
        out[f"dtol{fac}_t_horizon"] = c2.t_horizon
    # (c) manifold step eps ten times smaller (different nominal trajectories, same tube)
    c3 = mc.prepare_case(orb, n_phases=N_PHASES, eps=1e-7)
    out["eps_pos_raw"] = mc.run_levels(c3, SIG, zero, Z)["raw"]
    out["eps_t_horizon"] = c3.t_horizon
    # (d) a different random seed (how big is pure sampling noise?)
    out["seed_pos_raw"] = mc.run_levels(case, SIG, zero, mc.unit_draws(N_SAMPLES, SEED + 1))["raw"]
    save(name, case, **out)
    print(f"robustness z_far={z:.3f} done in {time.time()-t1:.0f}s", flush=True)

# ---- joint position x velocity grid for one orbit ---------------------------
name = "grid_z0.100.npz"
if not os.path.exists(os.path.join(args.out, name)) and not args.quick:
    t1 = time.time()
    orb = member_with_z(family, 0.10)
    case = mc.prepare_case(orb, n_phases=8)
    g = np.logspace(-9, -4, 11)
    SR, SV = np.meshgrid(g, g, indexing="ij")
    r = mc.run_levels(case, SR.ravel(), SV.ravel(), mc.unit_draws(200, SEED))
    np.savez_compressed(os.path.join(args.out, name), g=g, A=r["A"], B=r["B"])
    print(f"grid done in {time.time()-t1:.0f}s", flush=True)

json.dump({"seed": SEED, "n_phases": N_PHASES, "n_samples": N_SAMPLES, "z_far": Z_FAR,
           "mu": MU, "lunit_km": LUNIT_EM_KM, "total_seconds": time.time() - t0},
          open(os.path.join(args.out, "run_info.json"), "w"), indent=1)
print(f"ALL DONE in {time.time()-t0:.0f}s", flush=True)
