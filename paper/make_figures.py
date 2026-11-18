"""Regenerate every figure of the paper.

    python paper/make_figures.py          # figures from existing data/mc results
    python scripts/reproduce_all.py       # ONE command: experiment + analysis + figures + report + PDF

Figures 1-3 and 8a are computed here from the physics core in a few seconds.
Figures 4-7 and 8b read data/mc_summary.json (written by scripts/analyze_experiment.py).
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from cislab.constants import MU_EARTH_MOON as MU, LUNIT_EM_KM as KM, TUNIT_EM_S, C_ND_EM, R_MOON_KM
from cislab.cr3bp import lagrange_points, omega, jacobi
from cislab.cr3bp.core import grad_omega
from cislab.cr3bp.hill import critical_jacobi, hill_slice
from cislab.cr3bp.families import build_family, member_with_z, family_table
from cislab.cr3bp.manifolds import manifold, branch_towards_secondary, section_curve, convergence_check
from cislab.cr3bp import fastprop as fp

OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)

# Palette: validated for colour-vision deficiency (adjacent-pair CVD dE >= 12);
# line style is always used as a second cue, so no figure relies on colour alone.
BLUE, ORANGE, MAGENTA, GREEN, INK, GREY = "#2b7bba", "#d9541e", "#a8328e", "#3a9d5d", "#14304a", "#7d8b96"
SEQ = LinearSegmentedColormap.from_list("amp", ["#9ec5f4", "#256abf", "#0d366b"])   # one hue, light -> dark
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "text.color": INK, "axes.linewidth": 0.6,
    "grid.color": "#d5dde2", "grid.linewidth": 0.5, "lines.linewidth": 1.4, "legend.frameon": False,
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42,
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200)
    plt.close(fig)
    print("wrote", name)


def fig1_hill():
    """Lagrange points and Hill regions at three energies."""
    L, crit = lagrange_points(MU), critical_jacobi(MU)
    Cs = [crit["L1"] + 0.03, 0.5 * (crit["L1"] + crit["L2"]), crit["L2"] - 0.012]
    names = ["Closed: C above C(L1)", "L1 gate open", "L1 and L2 gates open"]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
    for ax, C, nm in zip(axs, Cs, names):
        u, v, W = hill_slice(MU, C, n=500)
        ax.contourf(u, v, W, levels=[-10, 0], colors=["#ead9ac"])
        ax.contour(u, v, W, levels=[0], colors=["#8a7335"], linewidths=0.7)
        ax.plot([-MU], [0], "o", ms=5, color=BLUE)
        ax.plot([1 - MU], [0], "o", ms=2.5, color=GREY)
        for k, p in L.items():
            ax.plot(p[0], p[1], "+", ms=6, mew=1, color=ORANGE if k in ("L1", "L2", "L3") else GREEN)
            ax.annotate(k, (p[0], p[1]), textcoords="offset points", xytext=(3, 4), fontsize=7)
        ax.set_title(f"{nm}\nC = {C:.4f}")
        ax.set_aspect("equal"); ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45); ax.set_xlabel("x")
    axs[0].set_ylabel("y")
    axs[0].annotate("Earth", (-MU, 0), textcoords="offset points", xytext=(4, -10), fontsize=7)
    axs[0].annotate("forbidden", (0.0, 1.02), fontsize=7, ha="center", color="#6b5a26")
    save(fig, "fig1_hill_regions")


def fig2_family(fam):
    tab = family_table(fam)
    z, T, k = np.abs(tab[:, 1]) * KM / 1e3, tab[:, 3] * TUNIT_EM_S / 86400, tab[:, 5]
    ipk = int(np.argmax(z))
    fig = plt.figure(figsize=(7.2, 3.0))
    ax = fig.add_subplot(1, 3, 1)
    pick = fam[4:ipk:14] + fam[ipk::28]
    for i, o in enumerate(pick):
        _, y = o.sample(300)
        ax.plot(y[0], y[2], color=SEQ(i / (len(pick) - 1)), lw=0.8)
    ax.plot([1 - MU], [0], "o", ms=3, color=GREY); ax.annotate("Moon", (1 - MU, 0), xytext=(3, -9), textcoords="offset points", fontsize=7)
    L2 = lagrange_points(MU)["L2"]
    ax.plot(L2[0], 0, "+", color=ORANGE, ms=6); ax.annotate("L2", (L2[0], 0), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax.set_xlabel("x"); ax.set_ylabel("z"); ax.set_aspect("equal"); ax.set_title("Family members, side view")
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.semilogy(z[:ipk], k[:ipk], color=BLUE, label="classic halo part")
    ax2.semilogy(z[ipk:], k[ipk:], color=ORANGE, ls="--", label="near-rectilinear part")
    ax2.axhline(1, color=GREY, lw=0.6); ax2.annotate("stable (k = 1)", (5, 1.08), fontsize=7, color=GREY)
    ax2.set_xlabel("z at far crossing (10$^3$ km)"); ax2.set_ylabel("stability index k"); ax2.grid(True); ax2.legend(fontsize=7)
    ax2.set_title("Instability falls along the family")
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(z[:ipk], T[:ipk], color=BLUE); ax3.plot(z[ipk:], T[ipk:], color=ORANGE, ls="--")
    # published catalogue values (NASA/JPL periodic-orbit API) for comparison
    jpl = [(0.12741504635570253, 3.2611173360604773), (0.026619595264470729, 3.4097239162419113), (0.18244450659922309, 1.5178035526683518)]
    ax3.plot([a * KM / 1e3 for a, _ in jpl], [b * TUNIT_EM_S / 86400 for _, b in jpl], "o", mfc="none", mec=INK, ms=6, label="JPL catalogue")
    ax3.set_xlabel("z at far crossing (10$^3$ km)"); ax3.set_ylabel("period (days)"); ax3.grid(True); ax3.legend(fontsize=7)
    ax3.set_title("Period, with published orbits")
    fig.tight_layout()
    save(fig, "fig2_halo_family")


def fig3_manifold(fam):
    o = member_with_z(fam, 0.10)
    sg = branch_towards_secondary(o, "stable")
    st = manifold(o, "stable", 60, 1e-6, sg, section_x=1 - MU)
    un = manifold(o, "unstable", 60, 1e-6, branch_towards_secondary(o, "unstable"), t_max=5.5)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9))
    ax = axs[0]
    for p in un:
        ax.plot(p.y[:, 0], p.y[:, 1], color=ORANGE, lw=0.35, alpha=0.7)
    for p in st:
        ax.plot(p.y[:, 0], p.y[:, 1], color=GREEN, lw=0.35, alpha=0.8)
    _, h = o.sample(300)
    ax.plot(h[0], h[1], color=INK, lw=1.4)
    ax.axvline(1 - MU, color=MAGENTA, lw=0.7, ls=":")
    ax.plot([], [], color=GREEN, label="stable (arrives)"); ax.plot([], [], color=ORANGE, label="unstable (leaves)")
    ax.legend(fontsize=6.5, loc="lower left"); ax.set_aspect("equal"); ax.set_xlim(0.9, 1.22); ax.set_ylim(-0.2, 0.2)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title("Moon-bound manifold branches")
    ph, cur = section_curve(st)
    _, cur10 = section_curve(manifold(o, "stable", 240, 1e-7, sg, section_x=1 - MU))
    for ax, (a, b, la, lb) in zip(axs[1:], [(1, 4, "y", "$v_y$"), (1, 2, "y", "z")]):
        ax.plot(cur10[:, a], cur10[:, b], color=GREY, lw=0.8, label=r"$\epsilon=10^{-7}$")
        ax.plot(cur[:, a], cur[:, b], "o", ms=2.6, color=MAGENTA, label=r"$\epsilon=10^{-6}$")
        ax.set_xlabel(la); ax.set_ylabel(lb); ax.grid(True)
    axs[1].legend(fontsize=6.5); axs[1].set_title("Capture curve at the Moon plane")
    th = np.linspace(0, 2 * np.pi, 100)
    axs[2].plot(R_MOON_KM / KM * np.cos(th), R_MOON_KM / KM * np.sin(th), color=INK, lw=0.8)
    axs[2].annotate("Moon's disc", (0, 0), ha="center", fontsize=6.5)
    axs[2].set_title("Same curve in position"); axs[2].set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    save(fig, "fig3_manifolds")
    return convergence_check(o, "stable", 1e-6, sg), o


def _amp_colors(halos):
    az = np.array([h["Az_km"] for h in halos])
    return [SEQ((a - az.min()) / max(az.max() - az.min(), 1)) for a in az]


def fig4_5_capture(S):
    H = S["halos"]
    cols = _amp_colors(H)
    for kind, xkey, xlabel, name, scale in (("pos", "sigma_km", "position error, 1$\\sigma$ per axis (m)", "fig4_capture_position", 1e3),
                                            ("vel", "sigma_mms", "velocity error, 1$\\sigma$ per axis (mm/s)", "fig5_capture_velocity", 1.0)):
        fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
        for ax, d, title in zip(axs, "AB", ["A: stays within 5000 km of the halo for one revolution",
                                             "B: stays within two Hill radii of the Moon"]):
            for h, c in zip(H, cols):
                ax.semilogx(np.array(h[xkey]) * scale, h[kind][d], color=c, lw=1.2)
            ax.axhline(0.5, color=GREY, lw=0.5, ls=":")
            ax.set_xlabel(xlabel); ax.set_title(title, fontsize=8); ax.grid(True)
        axs[0].set_ylabel("fraction captured")
        sm = plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(H[0]["Az_km"] / 1e3, H[-1]["Az_km"] / 1e3))
        cb = fig.colorbar(sm, ax=axs, pad=0.02, fraction=0.035)
        cb.set_label("halo amplitude $A_z$ (10$^3$ km)")
        save(fig, name)


def fig6_scaling(S):
    H = S["halos"]
    az = np.array([h["Az_km"] for h in H]) / 1e3
    k = np.array([h["stability"] for h in H])
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.8))
    for ax, kind, unit, scale, ylab in ((axs[0], "pos", "km", 1e3, "50% tolerance, position (m)"),
                                        (axs[1], "vel", "mms", 1.0, "50% tolerance, velocity (mm/s)")):
        for d, col, mk, ls in (("A", BLUE, "o", "-"), ("B", ORANGE, "s", "--")):
            v = np.array([h["s50"][kind][d][unit] for h in H]) * scale
            lo = np.array([h["s50"][kind][d]["lo"] for h in H]) * scale
            hi = np.array([h["s50"][kind][d]["hi"] for h in H]) * scale
            ax.errorbar(az, v, yerr=[v - lo, hi - v], color=col, marker=mk, ms=3.5, ls=ls, lw=1.1, capsize=1.5, elinewidth=0.6,
                        label=f"definition {d}")
        rs = os.path.join(ROOT, "data", "route_study.json")
        if os.path.exists(rs):          # linear prediction over ALL usable routes of a finer scan, 24 halos
            R = [r for r in json.load(open(rs))["linear"] if "pos_m" in r]
            key = "pos_m" if kind == "pos" else "vel_mms"
            ax.plot([r["Az_km"] / 1e3 for r in R], [r[key] for r in R], color=INK, lw=0.9, ls=":", label="linear theory, all routes (A)")
        else:
            lin = np.array([h["s50"][kind]["A_linear"] for h in H]) * scale
            ax.plot(az, lin, color=INK, lw=0.8, ls=":", label="linear theory (A)")
        ax.set_yscale("log"); ax.set_xlabel("halo amplitude $A_z$ (10$^3$ km)"); ax.set_ylabel(ylab); ax.grid(True)
    axs[0].legend(fontsize=6.5)
    ax = axs[2]
    v = np.array([h["s50"]["pos"]["A"]["km"] for h in H]) * 1e3
    ax.loglog(k, v, "o", color=BLUE, ms=3.5)
    f = S["fits"].get("pos_A")
    if f:
        kk = np.geomspace(k.min(), k.max(), 50)
        ax.loglog(kk, f["prefactor"] * kk ** f["exponent_vs_k"] * 1e3, color=INK, lw=0.8, ls=":")
        ax.annotate(f"slope {f['exponent_vs_k']:.2f}", (kk[len(kk) // 2], f["prefactor"] * kk[len(kk) // 2] ** f["exponent_vs_k"] * 1e3),
                    xytext=(6, 6), textcoords="offset points", fontsize=7)
    ax.set_xlabel("stability index k"); ax.set_ylabel("50% tolerance, position, A (m)"); ax.grid(True)
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout()
    save(fig, "fig6_tolerance_scaling")


def fig7_phase_and_grid(S):
    H = S["halos"]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.8))
    ax = axs[0]
    cols = _amp_colors(H)
    for h, c in list(zip(H, cols))[::2]:
        o = np.argsort(h["phases"])
        ax.semilogy(np.array(h["phases"])[o], np.array(h["s50"]["pos"]["A_by_phase"])[o] * 1e3, "o-", ms=2.5, lw=0.9, color=c,
                    label=f"{h['Az_km']/1e3:.0f}k km")
    ax.set_xlabel("arrival phase on the halo"); ax.set_ylabel("50% tolerance, position, A (m)")
    ax.legend(fontsize=6, title="$A_z$", title_fontsize=6.5); ax.grid(True); ax.set_title("Tolerance depends on the route")
    ax = axs[1]
    g = np.concatenate([h["gain_r"] for h in H]); s = np.concatenate([h["s50"]["pos"]["A_by_phase"] for h in H]) * 1e3
    ok = np.isfinite(s)
    ax.loglog(g[ok], s[ok], "o", ms=2.5, color=BLUE, alpha=0.8)
    gg = np.geomspace(g.min(), g.max(), 20)
    ax.loglog(gg, H[0]["d_tol_km"] * 1e3 / (0.6745 * gg), color=INK, ls=":", lw=0.8, label="$d_{tol}$ / (0.674 g)")
    ax.set_xlabel("error growth g (from the STM)"); ax.set_ylabel("50% tolerance per route (m)"); ax.legend(fontsize=6.5); ax.grid(True, which="both")
    ax.set_title("...and is set by error growth")
    ax = axs[2]
    if "grid" in S:
        G = S["grid"]
        x, y = np.array(G["sigma_km"]) * 1e3, np.array(G["sigma_mms"])
        pc = ax.pcolormesh(x, y, np.array(G["A"]).T, cmap=SEQ, vmin=0, vmax=1, shading="nearest")
        ax.contour(x, y, np.array(G["A"]).T, levels=[0.5], colors=[INK], linewidths=0.9)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("position error (m)"); ax.set_ylabel("velocity error (mm/s)")
        cb = fig.colorbar(pc, ax=ax, pad=0.02); cb.set_label("fraction captured (A)")
        ax.set_title("Joint errors, $A_z$ = 39k km")
    else:
        ax.set_axis_off()
    fig.tight_layout()
    save(fig, "fig7_routes_and_joint_errors")


def fig8_relativity(S, fam):
    o = member_with_z(fam, 0.10)
    t, y = o.sample(400)
    ratio = []
    for s in y.T:
        a_pn = np.array(fp.pn_test_particle(*s, MU, C_ND_EM))
        a_n = grad_omega(*s[:3], MU) - np.array([s[0], s[1], 0])          # gravity only (no centrifugal term)
        ratio.append(np.linalg.norm(a_pn) / np.linalg.norm(a_n))
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7))
    ax = axs[0]
    ax.semilogy(t / o.period, ratio, color=MAGENTA)
    ax.set_xlabel("position around the halo (revolutions)"); ax.set_ylabel("|1PN acceleration| / |Newtonian|"); ax.grid(True, which="both")
    ax.set_title("Size of the correction"); ax.set_ylim(1e-12, 1e-9)
    ax = axs[1]
    H = S["halos"]
    az = np.array([h["Az_km"] for h in H]) / 1e3
    eq = np.array([np.median(h["pn_equiv_error_m"]) for h in H])
    tol = np.array([np.nanmedian(h["s50"]["pos"]["A_by_phase"]) for h in H]) * 1e3
    ax.semilogy(az, tol, "o-", ms=3.5, color=BLUE, label="tolerance (median route)")
    ax.semilogy(az, eq, "s--", ms=3.5, color=MAGENTA, label="1PN, as an equivalent\ninjection error")
    ax.set_xlabel("halo amplitude $A_z$ (10$^3$ km)"); ax.set_ylabel("metres"); ax.legend(fontsize=6.5); ax.grid(True)
    ax.set_title("Relativity against the tolerance")
    ax = axs[2]
    h = H[len(H) // 2]
    x = np.array(h["sigma_km"]) * 1e3
    ax.semilogx(x, h["pos"]["A"], color=BLUE, lw=2.4, label="Newtonian")
    ax.semilogx(x, h["pos"]["A_pn"], color=MAGENTA, lw=0.9, ls="--", label="with 1PN")
    ax.set_xlabel("position error (m)"); ax.set_ylabel("fraction captured (A)"); ax.legend(fontsize=6.5); ax.grid(True)
    ax.set_title(f"{h['pos']['A_pn_flips']} of {h['pos']['A_n']} outcomes change")
    fig.tight_layout()
    save(fig, "fig8_relativity")
    return float(np.median(ratio)), float(np.max(ratio))


def main():
    fam = build_family(MU)
    fig1_hill()
    fig2_family(fam)
    conv, o = fig3_manifold(fam)
    extra = {"manifold_convergence": conv, "family_size": len(fam)}
    p = os.path.join(ROOT, "data", "mc_summary.json")
    if os.path.exists(p):
        S = json.load(open(p))
        fig4_5_capture(S); fig6_scaling(S); fig7_phase_and_grid(S)
        extra["pn_ratio_median"], extra["pn_ratio_max"] = fig8_relativity(S, fam)
    else:
        print("data/mc_summary.json missing: figures 4-8 skipped (run the experiment first)")
    json.dump(extra, open(os.path.join(ROOT, "data", "figure_numbers.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
