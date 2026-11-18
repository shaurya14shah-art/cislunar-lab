"""Turn the raw Monte Carlo tables (data/mc/*.npz) into the numbers the paper quotes.

    python scripts/analyze_experiment.py

Writes data/mc_summary.json.  Nothing here integrates anything: it only
classifies and counts, so it runs in seconds and can be re-run with different
capture definitions without repeating the expensive propagation.
"""
import glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cislab.cr3bp import montecarlo as mc
from cislab.constants import LUNIT_EM_KM, VUNIT_EM_KMS, TUNIT_EM_S

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "mc")
KM = LUNIT_EM_KM
MS = VUNIT_EM_KMS * 1000.0            # one velocity unit in m/s


def labels(raw, period, t_horizon, definition, n_rev=1.0):
    """raw: (levels, phases, samples, 7) -> labels (levels, phases, samples)."""
    L, P, S, _ = raw.shape
    out = np.empty((L, P, S), dtype=np.int8)
    for p in range(P):
        th = t_horizon[p] + (n_rev - 1.0) * period
        out[:, p] = mc.classify(raw[:, p].reshape(-1, 7).astype(float), period, th, n_rev, definition).reshape(L, S)
    return out


def curve(lab):
    return mc.capture_probability(lab)


def s50(sig, lab, unit):
    v = mc.sigma50(sig, curve(lab))
    lo, hi = mc.sigma50_bootstrap(sig, lab)
    return {"nd": v, "lo_nd": lo, "hi_nd": hi, unit[0]: v * unit[1], "lo": lo * unit[1], "hi": hi * unit[1]}


def summarize_halo(f):
    d = np.load(f)
    sig, T, th = d["sig"], float(d["period"]), d["t_horizon"]
    out = {"file": os.path.basename(f), "z_far": float(d["state0"][2]), "Az_km": float(d["amp"][2]) * KM,
           "Ax_km": float(d["amp"][0]) * KM, "Ay_km": float(d["amp"][1]) * KM,
           "period": T, "period_days": T * TUNIT_EM_S / 86400, "stability": float(d["stability"]),
           "jacobi": float(d["jacobi"]), "sigma_nd": sig.tolist(), "sigma_km": (sig * KM).tolist(),
           "sigma_mms": (sig * MS * 1000).tolist(), "phases": d["phases"].tolist(),
           "gain_r": d["gains"][:, 0].tolist(), "gain_v": d["gains"][:, 1].tolist(),
           "t_flight_days": (d["t_flight"] * TUNIT_EM_S / 86400).tolist(),
           "scan": d["scan"].tolist(), "d_tol_km": float(d["d_tol"]) * KM, "R_B_km": float(d["R_B"]) * KM,
           "pn_shift_m": (d["pn_shift"] * KM * 1000).tolist(),
           # the same shift expressed as the injection position error that would cause it
           "pn_equiv_error_m": (d["pn_shift"] / d["gains"][:, 0] * KM * 1000).tolist(),
           "s50": {}, "fractions": {}}
    for kind, unit in (("pos", ("km", KM)), ("vel", ("mms", MS * 1000))):
        raw, rawp = d[kind + "_raw"], d[kind + "_pn_raw"]
        out[kind] = {}
        out["s50"][kind] = {}
        for defn in "AB":
            lab, labp = labels(raw, T, th, defn), labels(rawp, T, th, defn)
            out[kind][defn] = curve(lab).tolist()
            out[kind][defn + "_pn"] = curve(labp).tolist()
            out[kind][defn + "_impact"] = (lab == mc.IMPACTED).mean(axis=(1, 2)).tolist()
            out[kind][defn + "_escape"] = (lab == mc.ESCAPED).mean(axis=(1, 2)).tolist()
            out["s50"][kind][defn] = s50(sig, lab, unit)
            out["s50"][kind][defn + "_pn"] = s50(sig, labp, unit)
            # per-phase tolerance: how much the answer depends on WHERE on the tube you inject
            out["s50"][kind][defn + "_by_phase"] = [mc.sigma50(sig, (lab[:, p] == mc.CAPTURED).mean(axis=1)) * unit[1]
                                                    for p in range(lab.shape[1])]
            # Newtonian vs 1PN with identical random draws: how many individual outcomes change?
            out[kind][defn + "_pn_flips"] = int((lab != labp).sum())
            out[kind][defn + "_n"] = int(lab.size)
            out[kind][defn + "_pn_maxdiff"] = float(np.max(np.abs(curve(lab) - curve(labp))))
        # sensitivity to the required dwell: 0.5, 1, 1.5, 2 revolutions (definition A)
        out["s50"][kind]["A_by_nrev"] = {str(n): mc.sigma50(sig, curve(labels(raw, T, th, "A", n))) * unit[1]
                                         for n in (0.5, 1.0, 1.5, 2.0)}
    # Linear-theory prediction (see montecarlo.prepare_case): a Gaussian error of size
    # sigma along the most-amplified direction is stretched by g; capture needs g*|xi| < d_tol,
    # which holds with probability 1/2 when sigma = d_tol / (0.6745 g).  The pooled curve is a
    # mixture over phases, so we predict it with the same mixture.
    for kind, col, unit in (("pos", 0, KM), ("vel", 1, MS * 1000)):
        g = d["gains"][:, col]
        from scipy.special import erf
        pred = np.array([np.mean(erf(float(d["d_tol"]) / (g * s * np.sqrt(2)))) for s in sig])
        out[kind]["A_linear"] = pred.tolist()
        out["s50"][kind]["A_linear"] = mc.sigma50(sig, pred) * unit
    return out


def summarize_robust(f):
    d = np.load(f)
    sig, T, th = d["sig"], float(d["period"]), d["t_horizon"]
    out = {"z_far": float(d["state0"][2]), "stability": float(d["stability"])}
    base = np.load(os.path.join(D, f"halo_z{out['z_far']:.3f}.npz"))
    for kind, unit in (("pos", KM),):
        for defn in "AB":
            a = labels(base[kind + "_raw"], T, base["t_horizon"], defn)
            b = labels(d[f"tight_{kind}_raw"], T, th, defn)
            out[f"tight_{kind}_{defn}"] = {
                "s50_base": mc.sigma50(sig, curve(a)) * unit, "s50_tight": mc.sigma50(sig, curve(b)) * unit,
                "max_curve_diff": float(np.max(np.abs(curve(a) - curve(b)))),
                "label_flips": int((a != b).sum()), "n": int(a.size)}
    a = labels(base["pos_raw"], T, base["t_horizon"], "A")
    out["s50_base_km"] = mc.sigma50(sig, curve(a)) * KM
    lo, hi = mc.sigma50_bootstrap(sig, a)
    out["s50_base_ci_km"] = [lo * KM, hi * KM]
    for tag, key in (("dtol0.5", "dtol0.5"), ("dtol2.0", "dtol2.0"), ("eps", "eps")):
        lab = labels(d[key + "_pos_raw"], T, d[key + "_t_horizon"], "A")
        out[f"{tag}_s50_km"] = mc.sigma50(sig, curve(lab)) * KM
    out["seed_s50_km"] = mc.sigma50(sig, curve(labels(d["seed_pos_raw"], T, th, "A"))) * KM
    return out


def main():
    files = sorted(glob.glob(os.path.join(D, "halo_z*.npz")))
    if not files:
        sys.exit("no results in data/mc - run scripts/run_experiment.py first")
    halos = [summarize_halo(f) for f in files]
    info = json.load(open(os.path.join(D, "run_info.json"))) if os.path.exists(os.path.join(D, "run_info.json")) else {}
    d0 = np.load(files[0])
    summary = {"halos": halos, "n_phases": int(d0["pos_raw"].shape[1]), "n_samples": int(d0["pos_raw"].shape[2]),
               "n_levels": int(d0["pos_raw"].shape[0]), "run_info": info, "robust": [], "fits": {}}
    for f in sorted(glob.glob(os.path.join(D, "robust_z*.npz"))):
        summary["robust"].append(summarize_robust(f))
    # Scaling law: fit  log(sigma50) = a + b log(k)  over the unstable members.
    k = np.array([h["stability"] for h in halos])
    az = np.array([h["Az_km"] for h in halos])
    for kind, unit in (("pos", "km"), ("vel", "mms")):
        for defn in "AB":
            s = np.array([h["s50"][kind][defn][unit] for h in halos])
            ok = np.isfinite(s) & (k > 1.5)
            if ok.sum() >= 3:
                b, a = np.polyfit(np.log(k[ok]), np.log(s[ok]), 1)
                res = np.log(s[ok]) - (a + b * np.log(k[ok]))
                b2, a2 = np.polyfit(np.log(az[ok]), np.log(s[ok]), 1)
                summary["fits"][f"{kind}_{defn}"] = {"exponent_vs_k": float(b), "prefactor": float(np.exp(a)),
                                                     "rms_log_residual": float(np.sqrt(np.mean(res ** 2))),
                                                     "exponent_vs_Az": float(b2), "n": int(ok.sum())}
    g = os.path.join(D, "grid_z0.100.npz")
    if os.path.exists(g):
        gd = np.load(g)
        n = len(gd["g"])
        summary["grid"] = {"sigma_km": (gd["g"] * KM).tolist(), "sigma_mms": (gd["g"] * MS * 1000).tolist(),
                           "A": mc.capture_probability(gd["A"]).reshape(n, n).tolist(),
                           "B": mc.capture_probability(gd["B"]).reshape(n, n).tolist()}
    json.dump(summary, open(os.path.join(ROOT, "data", "mc_summary.json"), "w"), indent=1)
    print(f"{'z_far':>6} {'Az km':>8} {'k':>7} | {'s50 pos A (m)':>14} {'B (m)':>10} | {'s50 vel A (mm/s)':>16} {'B':>9} | PN flips A")
    for h in halos:
        print(f"{h['z_far']:6.3f} {h['Az_km']:8.0f} {h['stability']:7.1f} | {h['s50']['pos']['A']['km']*1e3:14.2f} "
              f"{h['s50']['pos']['B']['km']*1e3:10.2f} | {h['s50']['vel']['A']['mms']:16.4f} {h['s50']['vel']['B']['mms']:9.4f} | "
              f"{h['pos']['A_pn_flips']}/{h['pos']['A_n']}")
    print("fits:", json.dumps(summary["fits"], indent=1))
    for r in summary["robust"]:
        print("robust", json.dumps(r))


if __name__ == "__main__":
    main()
