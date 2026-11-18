"""Run the ten mandatory checks, print measured values, and write docs/VERIFICATION.md.

    python scripts/verification_report.py

Every line of the report is produced by executing code now - nothing is copied
from an earlier run.  The pass/fail thresholds are the same as in
tests/test_verification.py.
"""
import datetime, json, os, platform, subprocess, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from cislab.constants import MU_EARTH_MOON as MU, GM_SUN, C, ARCSEC
from cislab.cr3bp import propagate, lagrange_points, linear_frequencies, RTOL, ATOL
from cislab.cr3bp.core import jacobi_drift
from cislab.cr3bp.families import build_family, member_with_z
from cislab.cr3bp import montecarlo as mc
from cislab.cr3bp.periodic import correct_fixed
from cislab.relativity import schwarzschild as sw
from cislab.relativity.clocks import gps_clock_offset
from cislab.relativity.eih import mercury_precession_arcsec_per_century

rows = []


def add(n, what, measured, reference, ok):
    rows.append((n, what, measured, reference, "PASS" if ok else "FAIL"))
    print(f"[{'PASS' if ok else 'FAIL'}] {n}. {what}: {measured}  (reference: {reference})")


# 1 ---------------------------------------------------------------------------
r0 = 0.7
s0 = np.array([r0, 0, 0, 0, r0 * (r0 ** -1.5 - 1), 0])
sol = propagate(s0, (0, 50), 0.0, t_eval=np.linspace(0, 50, 2001))
dev = float(np.max(np.abs(np.linalg.norm(sol.y[:3], axis=0) - r0)))
add(1, "mu = 0: circular orbit stays circular over 50 time units (max radius change)", f"{dev:.1e}", "< 1e-12", dev < 1e-12)
# 2 ---------------------------------------------------------------------------
L = lagrange_points(MU)
sol = propagate(np.array([*L["L4"], 0, 0, 0]), (0, 100), MU)
dev = float(np.max(np.abs(sol.y[:3].T - L["L4"])))
add(2, "Particle at L4 with zero velocity, 100 time units (max displacement)", f"{dev:.1e}", "< 1e-11", dev < 1e-11)
# 3 ---------------------------------------------------------------------------
ref = {"L1": 0.836915, "L2": 1.155682, "L3": -1.005063}
ok = all(abs(L[k][0] - v) < 5e-7 for k, v in ref.items())
add(3, "Collinear points x(L1), x(L2), x(L3)", ", ".join(f"{L[k][0]:.9f}" for k in ref), "0.836915, 1.155682, -1.005063", ok)
# 4 ---------------------------------------------------------------------------
f1, f2 = linear_frequencies(MU, "L1"), linear_frequencies(MU, "L2")
ok = np.allclose(f1, (2.9321, 2.3344, 2.2688), atol=6e-5) and np.allclose(f2, (2.1587, 1.8626, 1.7862), atol=6e-5)
add(4, "Linear rates (lambda, nu, omega_z) at L1 and L2", "L1 " + ", ".join(f"{v:.5f}" for v in f1) + "; L2 " + ", ".join(f"{v:.5f}" for v in f2),
    "L1 2.9321, 2.3344, 2.2688; L2 2.1587, 1.8626, 1.7862", ok)
t = np.linspace(2.0, 4.5, 60)
sol = propagate(np.array([L["L1"][0] + 1e-8, 0, 0, 0, 0, 0]), (0, 4.5), MU, t_eval=t)
slope = float(np.polyfit(t, np.log(np.linalg.norm(sol.y[:3].T - L["L1"], axis=1)), 1)[0])
add(4, "Measured divergence rate of a 1e-8 displacement at L1", f"{slope:.5f}", "lambda = 2.9321 (within 0.2%)", abs(slope - 2.9321) / 2.9321 < 2e-3)
# 5 ---------------------------------------------------------------------------
fam = build_family(MU)
orb = member_with_z(fam, 0.10)
d1 = jacobi_drift(propagate(orb.state0, (0, 5 * orb.period), MU), MU)
d2 = jacobi_drift(propagate(np.array([0.3, 0.0, 0.05, 0.0, 1.9, 0.0]), (0, 60), MU), MU)
add(5, "Relative Jacobi drift: 5 halo revolutions / 60 time-unit free trajectory", f"{d1:.1e} / {d2:.1e}", "< 1e-11", max(d1, d2) < 1e-11)
# 6 ---------------------------------------------------------------------------
JPL = [(1.1589028132163817, 0.12741504635570253, -0.20997144894179126, 3.0901151159191, 3.2611173360604773, 170.405268543231, "first"),
       (1.1801895376953091, 0.026619595264470729, -0.15983377173327767, 3.14903225455325, 3.4097239162419113, 574.471019364258, "first"),
       (1.0225293572044083, 0.18244450659922309, -0.10435995011600951, 3.04599801919558, 1.5178035526683518, 1.33678978716921, "last")]
for x0, z0, vy0, Cj, T, stab, br in JPL:
    o = member_with_z(fam, z0, br)
    errs = [abs(o.state0[0] - x0) / x0, abs(o.state0[4] - vy0) / abs(vy0), abs(o.period - T) / T,
            abs(o.jacobi - MU * (1 - MU) - Cj) / Cj]
    add(6, f"JPL catalogue L2 halo with z0 = {z0:.6f}: x0, vy0, T, C_JPL",
        f"{o.state0[0]:.10f}, {o.state0[4]:.10f}, {o.period:.10f}, {o.jacobi - MU*(1-MU):.10f} (k = {o.stability_index:.4f})",
        f"{x0:.10f}, {vy0:.10f}, {T:.10f}, {Cj:.10f} (k = {stab:.4f}); worst relative error {max(errs):.1e}", max(errs) < 1e-7)
# 7 ---------------------------------------------------------------------------
m = mercury_precession_arcsec_per_century(True)
m0 = mercury_precession_arcsec_per_century(False)
add(7, "Mercury perihelion advance from EIH integration (control with 1PN off)", f"{m:.4f} arcsec/century ({m0:.4f})", "42.98 within 1%", abs(m - 42.98) / 42.98 < 0.01)
# 8 ---------------------------------------------------------------------------
isco = sw.isco_numeric()
add(8, "Schwarzschild ISCO found numerically / photon sphere", f"{isco:.9f} / {sw.photon_sphere_numeric():.9f} GM/c^2", "6 / 3", abs(isco - 6) < 1e-6)
b = 6.957e8 / (GM_SUN / C ** 2)
dfl = sw.light_deflection(b) / ARCSEC
add(8, "Light deflection at the solar limb: exact null geodesic (weak-field formula)", f"{dfl:.5f} arcsec ({sw.light_deflection_weak(b)/ARCSEC:.5f})", "about 1.75", abs(dfl - 1.75) < 0.005)
# 9 ---------------------------------------------------------------------------
g = gps_clock_offset()
add(9, "GPS net clock offset", f"+{g['microseconds_per_day']:.3f} microseconds/day (fraction {g['fractional_offset']:.4e})",
    "about +38; Ashby (2003) 4.4647e-10", abs(g["microseconds_per_day"] - 38) < 1)
# 10 --------------------------------------------------------------------------
X = np.array([orb.state0[0], orb.state0[2], orb.state0[4], orb.period / 2])
o2, _ = correct_fixed(X, MU, "z0", rtol=RTOL / 2, atol=ATOL / 2)
dT = abs(o2.period - orb.period)
case = mc.prepare_case(orb, n_phases=4)
Z = mc.unit_draws(120)
sig = np.logspace(-9, -4, 11)
a = mc.run_levels(case, sig, 0 * sig, Z)
bb = mc.run_levels(case, sig, 0 * sig, Z, rtol=5e-13, atol=5e-15, hmax=0.005)
dp = max(float(np.max(np.abs(mc.capture_probability(a[d]) - mc.capture_probability(bb[d])))) for d in "AB")
flips = int((a["A"] != bb["A"]).sum())
add(10, "Halving every tolerance: halo period change / largest change in any capture probability / outcomes changed",
    f"{dT:.1e} / {dp:.4f} / {flips} of {a['A'].size}", "no conclusion changes (< 1e-11, < 0.03)", dT < 1e-11 and dp < 0.03)
summ = os.path.join(ROOT, "data", "mc_summary.json")
if os.path.exists(summ):
    for r in json.load(open(summ)).get("robust", []):
        t_ = r["tight_pos_A"]
        add(10, f"Full batch run, halo z = {r['z_far']:.2f}: 50% position tolerance, normal vs halved tolerances",
            f"{t_['s50_base']*1e3:.2f} m vs {t_['s50_tight']*1e3:.2f} m ({t_['label_flips']} of {t_['n']} outcomes changed)",
            "same within sampling interval", abs(np.log(t_["s50_base"] / t_["s50_tight"])) < 0.1)

# --- full test suite -----------------------------------------------------------
res = subprocess.run([sys.executable, "-m", "pytest", "-q", os.path.join(ROOT, "tests")], capture_output=True, text=True, cwd=ROOT)
tail = res.stdout.strip().splitlines()[-1] if res.stdout.strip() else res.stderr.strip()[-200:]
print("pytest:", tail)

import numpy, scipy
try:
    import numba; nb = numba.__version__
except Exception:
    nb = "not installed"
with open(os.path.join(ROOT, "docs", "VERIFICATION.md"), "w") as f:
    f.write("# Verification report\n\n")
    f.write(f"Generated by `scripts/verification_report.py` on {datetime.date.today().isoformat()} "
            f"(Python {platform.python_version()}, NumPy {numpy.__version__}, SciPy {scipy.__version__}, Numba {nb}).\n\n")
    f.write("Every number in the *Measured* column was computed when this file was generated. "
            "The numbering follows section 6 of the project brief.\n\n")
    f.write("| # | Check | Measured | Reference / threshold | Result |\n|---|---|---|---|---|\n")
    for n, what, meas, ref_, ok in rows:
        f.write(f"| {n} | {what} | {meas} | {ref_} | **{ok}** |\n")
    f.write(f"\nFull automated test suite (`pytest tests`): **{tail}**\n\n")
    f.write("## What these checks do not cover\n\n"
            "* The 1PN *energy* formula was not read from a primary source; it is validated only by the conservation test "
            "`test_1pn_energy_is_conserved_better_than_newtonian_energy`.\n"
            "* Richardson's coefficients were cross-checked against an open-source implementation, not against the 1980 article "
            "(paywalled). They are validated indirectly: the seed converges to orbits that match JPL's catalogue.\n"
            "* The browser interface was exercised with automated screenshots of the main interactions, not with a full "
            "interaction test suite.\n"
            "* Nothing here validates the CR3BP *as a model of the real Earth-Moon system*; see the paper's limitations section.\n")
print("wrote docs/VERIFICATION.md;", sum(r[-1] == "PASS" for r in rows), "of", len(rows), "checks pass")
sys.exit(0 if all(r[-1] == "PASS" for r in rows) and res.returncode == 0 else 1)
