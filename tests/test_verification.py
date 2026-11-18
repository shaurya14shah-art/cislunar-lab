"""The ten mandatory verification checks (brief, section 6), numbered as in the brief.

Run:  pytest tests/test_verification.py -v
Each test states what is being checked, against which reference, and how tight
the agreement must be.  docs/VERIFICATION.md is generated from the same checks
by scripts/verification_report.py.
"""
import numpy as np
import pytest
from cislab.cr3bp import (propagate, lagrange_points, linear_frequencies, jacobi, RTOL, ATOL)
from cislab.cr3bp.core import jacobi_drift, rhs
from cislab.cr3bp.families import member_with_z
from cislab.cr3bp import fastprop as fp
from cislab.constants import ARCSEC, GM_SUN, C

# --- published references ---------------------------------------------------
# NASA/JPL Three-Body Periodic Orbits API v1.0, Earth-Moon, L2 halo, branch N
# (queries made 10 Oct 2026; see docs/SOURCES.md).  Columns: x0, z0, vy0, C, T, stability
JPL_HALOS = [
    (1.1589028132163817, 0.12741504635570253, -0.20997144894179126, 3.0901151159191,
     3.2611173360604773, 170.405268543231),
    (1.1801895376953091, 0.026619595264470729, -0.15983377173327767, 3.14903225455325,
     3.4097239162419113, 574.471019364258),
    (1.0225293572044083, 0.18244450659922309, -0.10435995011600951, 3.04599801919558,
     1.5178035526683518, 1.33678978716921),
]


def test_01_mu_zero_is_kepler():
    """mu = 0: one primary at the origin.  (a) A circular orbit stays circular.
    (b) An eccentric orbit conserves inertial energy and angular momentum and
    returns to its start after one Kepler period 2*pi*a^1.5."""
    r0 = 0.7
    n = r0 ** -1.5                                  # inertial angular speed of a circular orbit
    s0 = np.array([r0, 0, 0, 0, r0 * (n - 1), 0])   # speed seen in the frame rotating at rate 1
    sol = propagate(s0, (0, 50), 0.0, t_eval=np.linspace(0, 50, 2001))
    r = np.linalg.norm(sol.y[:3], axis=0)
    assert np.max(np.abs(r - r0)) < 1e-12
    # (b) eccentric: start at periapsis rp with e = 0.4
    e, a = 0.4, 0.8
    rp = a * (1 - e)
    vp = np.sqrt((1 + e) / rp)
    T = 2 * np.pi * a ** 1.5
    s0 = np.array([rp, 0, 0, 0, vp - rp, 0])
    sf = propagate(s0, (0, T), 0.0).y[:, -1]

    def inertial(s, t):
        c, s_ = np.cos(t), np.sin(t)
        R = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
        return R @ s[:3], R @ (s[3:] + np.cross([0, 0, 1], s[:3]))
    p0, v0 = inertial(s0, 0.0)
    p1, v1 = inertial(sf, T)
    assert np.linalg.norm(p1 - p0) < 1e-11 and np.linalg.norm(v1 - v0) < 1e-11


def test_02_L4_is_an_equilibrium(mu):
    """A particle at rest at L4 stays there for 100 time units (~16 months)."""
    L4 = lagrange_points(mu)["L4"]
    sol = propagate(np.array([*L4, 0, 0, 0]), (0, 100), mu)
    assert np.max(np.abs(sol.y[:3].T - L4)) < 1e-11
    assert np.max(np.abs(sol.y[3:])) < 1e-11


def test_03_collinear_points_six_figures(mu):
    L = lagrange_points(mu)
    for name, ref in (("L1", 0.836915), ("L2", 1.155682), ("L3", -1.005063)):
        assert abs(L[name][0] - ref) < 5e-7, name
    # and against JPL's 15-digit values
    assert abs(L["L1"][0] - 0.836915125772357) < 1e-13
    assert abs(L["L2"][0] - 1.15568216544488) < 1e-13


def test_04_linear_rates_and_measured_divergence(mu):
    ref = {"L1": (2.9321, 2.3344, 2.2688), "L2": (2.1587, 1.8626, 1.7862)}
    for p, vals in ref.items():
        got = linear_frequencies(mu, p)
        assert np.allclose(got, vals, atol=6e-5), (p, got)
    # Measured growth: displace by 1e-8 at L1 and fit the slope of ln(distance).
    L1 = lagrange_points(mu)["L1"]
    s0 = np.array([L1[0] + 1e-8, 0, 0, 0, 0, 0])
    t = np.linspace(2.0, 4.5, 60)          # after the transient, before nonlinearity (d < 1e-2)
    sol = propagate(s0, (0, 4.5), mu, t_eval=t)
    d = np.linalg.norm(sol.y[:3].T - L1, axis=1)
    assert d.max() < 1e-2
    slope = np.polyfit(t, np.log(d), 1)[0]
    assert abs(slope - 2.9321) / 2.9321 < 2e-3, slope


def test_05_jacobi_drift(mu, family):
    """Relative Jacobi drift < 1e-11 over five halo revolutions and over a
    60-time-unit (~9.5 revolutions of the Moon) transit trajectory."""
    orb = member_with_z(family, 0.10)
    sol = propagate(orb.state0, (0, 5 * orb.period), mu)
    assert jacobi_drift(sol, mu) < 1e-11
    sol = propagate(np.array([0.3, 0.0, 0.05, 0.0, 1.9, 0.0]), (0, 60), mu)
    assert jacobi_drift(sol, mu) < 1e-11


@pytest.mark.parametrize("ref", JPL_HALOS)
def test_06_published_halo(mu, family, ref):
    """Reproduce JPL catalogue halos starting from OUR family (Richardson seed ->
    differential correction -> continuation); only z0 is taken from the catalogue.
    x0, vy0, period and Jacobi constant must then agree to >= 7 significant figures
    (the brief asks for 4)."""
    x0, z0, vy0, Cj, T, stab = ref
    branch = "last" if T < 2.0 else "first"
    orb = member_with_z(family, z0, branch)
    assert abs(orb.state0[0] - x0) / x0 < 1e-7
    assert abs(orb.state0[4] - vy0) / abs(vy0) < 1e-7
    assert abs(orb.period - T) / T < 1e-7
    # JPL's Jacobi constant omits the constant mu(1-mu); convert before comparing.
    assert abs(orb.jacobi - mu * (1 - mu) - Cj) / Cj < 1e-9
    assert abs(orb.stability_index - stab) / stab < 1e-4      # same stability definition as JPL


def test_07_mercury_perihelion():
    from cislab.relativity.eih import mercury_precession_arcsec_per_century
    got = mercury_precession_arcsec_per_century(pn=True)
    assert abs(got - 42.98) / 42.98 < 0.01, got
    # control: the same code with relativity switched off gives no precession
    assert abs(mercury_precession_arcsec_per_century(pn=False)) < 0.01


def test_08_schwarzschild_isco_and_light_bending():
    from cislab.relativity import schwarzschild as sw
    assert abs(sw.isco_numeric() - 6.0) < 1e-6
    assert abs(sw.photon_sphere_numeric() - 3.0) < 1e-6
    # dynamical check of the ISCO: a slightly perturbed circular orbit survives
    # at r = 6.5 M but plunges from r = 5.5 M
    for r, plunges in ((6.5, False), (5.5, True)):
        L, E = np.sqrt(sw.circular_L2(r)), sw.circular_E(r)
        g = sw.geodesic(r, 0.0, E * (1 + 1e-9), L, 3000.0, inward=True)
        assert g["plunged"] == plunges
    b = 6.957e8 / (GM_SUN / C ** 2)          # solar radius in units of GM/c^2
    defl = sw.light_deflection(b) / ARCSEC
    assert abs(defl - 1.75) < 0.005, defl
    assert abs(sw.light_deflection_weak(b) / ARCSEC - defl) < 1e-4   # weak-field formula check


def test_09_gps_clock():
    from cislab.relativity.clocks import gps_clock_offset
    r = gps_clock_offset()
    assert abs(r["microseconds_per_day"] - 38.0) < 1.0
    assert abs(r["fractional_offset"] - 4.4647e-10) < 5e-14      # Ashby (2003) Eq. 35


def test_10_halving_tolerances_changes_nothing(mu, family):
    """(a) the corrected halo is the same at half the tolerance;
       (b) a small Monte Carlo gives the same capture curve and the same
           50%-tolerance when rtol, atol and the step cap are all halved."""
    from cislab.cr3bp.periodic import correct_fixed
    from cislab.cr3bp import montecarlo as mc
    orb = member_with_z(family, 0.10)
    X = np.array([orb.state0[0], orb.state0[2], orb.state0[4], orb.period / 2])
    X[[0, 2, 3]] *= 1 + 1e-6                     # start Newton slightly off; z0 is held fixed
    o2, _ = correct_fixed(X, mu, "z0", rtol=RTOL / 2, atol=ATOL / 2)
    assert abs(o2.period - orb.period) < 1e-11 and np.max(np.abs(o2.state0 - orb.state0)) < 1e-11
    assert abs(o2.stability_index - orb.stability_index) / orb.stability_index < 1e-7

    case = mc.prepare_case(orb, n_phases=4)
    Z = mc.unit_draws(120)
    sig = np.logspace(-9, -4, 11)
    zero = np.zeros_like(sig)
    a = mc.run_levels(case, sig, zero, Z)
    b = mc.run_levels(case, sig, zero, Z, rtol=5e-13, atol=5e-15, hmax=0.005)
    for d in "AB":
        pa, pb = mc.capture_probability(a[d]), mc.capture_probability(b[d])
        assert np.max(np.abs(pa - pb)) < 0.03, (d, pa, pb)
    sa, sb = mc.sigma50(sig, mc.capture_probability(a["A"])), mc.sigma50(sig, mc.capture_probability(b["A"]))
    assert abs(np.log(sa / sb)) < 0.1
