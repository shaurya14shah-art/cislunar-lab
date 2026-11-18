"""Tests for the N-body sandbox, the relativity layer and the field tools."""
import numpy as np
import pytest
from cislab.constants import (PRESETS, GM_EARTH, GM_MOON, GM_SUN, C, G, Body, MU_EARTH_MOON, C_ND_EM)
from cislab.nbody.system import NBodySystem, two_body
from cislab.nbody.kernels import accel_newton, accel_eih_correction
from cislab.relativity import regime, schwarzschild as sw, clocks, light, eih
from cislab import fields
from cislab.cr3bp import fastprop as fp


def test_two_body_period_and_symplectic_energy():
    s = two_body()
    P = 2 * np.pi * np.sqrt(3.844e8 ** 3 / (GM_EARTH + GM_MOON))
    h = s.run(50 * P, dt=P / 1000, method="yoshida4", update=False)
    assert np.max(np.abs(h["energy_drift"])) < 1e-10          # bounded, no secular drift
    h1 = s.run(P, dt=P / 4000, method="yoshida4", update=False)
    assert np.max(np.abs(h1["pos"][-1] - h1["pos"][0])) / 3.844e8 < 1e-8   # closes after one period
    assert not h["collisions"]


def test_yoshida_is_fourth_order():
    s = two_body()
    P = 2 * np.pi * np.sqrt(3.844e8 ** 3 / (GM_EARTH + GM_MOON))
    ref = s.run(P, method="dop853", rtol=1e-13, n_frames=1, update=False)["pos"][-1]
    err = [np.max(np.abs(s.run(P, dt=P / n, method="yoshida4", n_frames=1, update=False)["pos"][-1] - ref))
           for n in (200, 400)]
    order = np.log2(err[0] / err[1])
    assert 3.6 < order < 4.4, order


def test_momentum_and_angular_momentum_conserved():
    s = NBodySystem([PRESETS["Sun"]])
    s.add_circular_satellite("Earth", 0, 1.496e11)
    s.add_circular_satellite("Jupiter", 0, 7.78e11)
    s.to_barycentric()
    h = s.run(3.2e8, method="dop853", update=False)
    p = np.array([s.momentum(v) for v in h["vel"]])
    Lz = np.array([s.angular_momentum(a, b)[2] for a, b in zip(h["pos"], h["vel"])])
    assert np.max(np.abs(p)) / (GM_SUN * 3e4) < 1e-12
    assert np.max(np.abs(Lz - Lz[0])) / abs(Lz[0]) < 1e-10


def test_collision_is_flagged_not_merged():
    s = NBodySystem([Body("A", GM_EARTH, 6.4e6, (0, 0, 0), (0, 0, 0)),
                     Body("B", GM_MOON, 1.7e6, (5e7, 0, 0), (0, 0, 0))])   # dropped from rest
    h = s.run(2e4, dt=5.0, method="yoshida4", update=False)
    assert len(h["collisions"]) == 1 and h["collisions"][0]["names"] == ["A", "B"]
    assert len(s.bodies) == 2                                    # nothing was merged


def test_corotating_frame_freezes_the_pair():
    s = two_body()
    h = s.run(2e6, method="dop853", update=False)
    Pr, Vr = s.transform(h, "corotating", (0, 1))
    assert np.max(np.abs(Pr[:, :, 1:])) < 1e-3 and np.max(np.abs(Vr)) < 1e-6
    Pb, Vb = s.transform(h, "barycentric")
    assert np.max(np.abs((s.gm[None, :, None] * Pb).sum(1))) / (GM_EARTH * 3.8e8) < 1e-12


def test_editing_and_units():
    s = two_body()
    k = s.add("Jupiter", pos=(1e12, 0, 0), vel=(0, 1e3, 0), name="J")
    s.edit(k, radius=5e7, mass=1e27)
    assert abs(s.bodies[k].gm - 1e27 * G) < 1 and s.bodies[k].radius == 5e7
    s.remove(k)
    assert len(s.bodies) == 2
    u = s.unit_system((0, 1))
    assert abs(u["mu"] - GM_MOON / (GM_EARTH + GM_MOON)) < 1e-15
    assert abs(u["time_s"] * 2 * np.pi / 86400 - 27.28) < 0.05   # sidereal month from Kepler's law
    with pytest.raises(ValueError):
        s.run(1e3, dt=1.0, method="yoshida4", pn=True)


def test_eih_reduces_to_newton_and_scales_as_c_minus_2():
    s = two_body("Sun", "Mercury", 5.79e10)
    aN = s.acceleration()
    a1 = s.acceleration(pn=True)
    s2 = NBodySystem(s.bodies, c=10 * C)
    a2 = s2.acceleration(pn=True)
    r1, r2 = np.abs(a1 - aN).max(), np.abs(a2 - aN).max()
    assert 1e-9 < r1 / np.abs(aN).max() < 1e-6
    assert abs(r1 / r2 - 100) < 1e-3                # correction is exactly proportional to 1/c^2


def test_1pn_energy_is_conserved_better_than_newtonian_energy():
    """A compact binary with v/c ~ 0.02.  Under EIH dynamics the *Newtonian*
    energy wobbles at relative size ~(v/c)^2 while the 1PN energy holds ~(v/c)^4.
    This validates both the force law and the energy formula against each other."""
    ns = PRESETS["Neutron star"]
    s = NBodySystem([Body("a", ns.gm, ns.radius), Body("b", 0.7 * ns.gm, ns.radius)])
    d = 2.0e7
    s.bodies[1].pos = (d, 0, 0)
    v = 0.85 * np.sqrt((s.bodies[0].gm + s.bodies[1].gm) / d)      # eccentric
    s.bodies[1].vel = (0, v, 0)
    s.to_barycentric()
    assert 1e-4 < (v / C) ** 2 < 1e-3
    P = 2 * np.pi * np.sqrt(d ** 3 / (s.gm.sum()))
    h = s.run(1.5 * P, method="dop853", pn=True, rtol=1e-13, update=False)
    eN = np.array([s.energy(p, w, pn=False) for p, w in zip(h["pos"], h["vel"])])
    wobN = np.ptp(eN) / abs(eN[0])
    wobP = np.max(np.abs(h["energy_drift"]))
    assert wobN > 1e-5 and wobP < wobN * 5e-3, (wobN, wobP)


def test_test_particle_1pn_matches_general_eih():
    """The specialised CR3BP 1PN term (rotating frame) must equal the general
    EIH kernel evaluated for Earth, Moon and a massless spacecraft."""
    mu, c = MU_EARTH_MOON, C_ND_EM
    x, y, z, vx, vy, vz = 1.08, 0.03, 0.12, 0.05, -0.2, 0.07
    got = np.array(fp.pn_test_particle(x, y, z, vx, vy, vz, mu, c))
    pos = np.array([[-mu, 0, 0], [1 - mu, 0, 0], [x, y, z]])
    vel = np.array([[0, -mu, 0], [0, 1 - mu, 0], [vx - y, vy + x, vz]])     # inertial velocities
    gm = np.array([1 - mu, mu, 0.0])
    aN = np.empty((3, 3)); accel_newton(pos, gm, aN)
    da = np.empty((3, 3)); accel_eih_correction(pos, vel, gm, c, aN, da)
    assert np.allclose(got, da[2], rtol=1e-12, atol=0)
    # brief's estimate: the correction is ~1e-11 of the Newtonian acceleration
    ratio = np.linalg.norm(got) / np.linalg.norm(aN[2])
    assert 1e-12 < ratio < 1e-9, ratio


def test_regime_meter():
    s = two_body()
    rep = regime.regime_report(s.gm, s.pos, s.vel, s.radius, ["Earth", "Moon"])
    moon = rep[1]
    assert 1e-12 < moon["eps_orbit"] < 2e-11            # GM_E/(r c^2) ~ 1.15e-11
    assert moon["band"] == "newtonian" and not moon["surface_warning"]
    ns = PRESETS["Neutron star"]
    rep = regime.regime_report([ns.gm, ns.gm], [[0, 0, 0], [1e5, 0, 0]], [[0, 0, 0], [0, 4e7, 0]],
                               [ns.radius, ns.radius])
    assert rep[0]["band"] in ("1pn-marginal", "breakdown", "1pn-ok") and rep[0]["surface_warning"]
    assert abs(regime.lorentz_gamma(0.6 * C) - 1.25) < 1e-12
    assert abs(regime.kinetic_energy_ratio(1.0) - 1) < 1e-12
    assert abs(regime.kinetic_energy_ratio(0.6 * C) - 0.25 / 0.18) < 1e-12


def test_schwarzschild_orbits_and_clocks():
    # weak-field limit of the exact precession matches 6 pi M / (a (1-e^2))
    a, e = 1e6, 0.3
    exact = sw.precession_per_orbit(a * (1 - e), a * (1 + e))
    assert abs(exact / (6 * np.pi / (a * (1 - e * e))) - 1) < 1e-4
    g = sw.bound_orbit(10.0, 30.0, n_orbits=3)
    assert not g["plunged"] and g["r"].min() > 9.99 and g["r"].max() < 30.01
    # proper time runs slower than far-away time
    assert g["t"][-1] > g["tau"][-1]
    assert abs(sw.clock_rate_static(4.0) - np.sqrt(0.5)) < 1e-15
    assert np.isnan(sw.light_deflection(5.0))            # below 3*sqrt(3): captured
    assert sw.light_ray(-40.0, 4.0, 0.0)["captured"] and not sw.light_ray(-40.0, 12.0, 0.0)["captured"]


def test_weak_field_ray_tracer_matches_formula():
    b = 7e8 * 5
    out = light.trace_ray([-4e12, b, 0], [1, 0, 0], [GM_SUN], [[0, 0, 0]], 8e12, n=20000)
    angle = np.arctan2(-out["direction"][1], out["direction"][0])
    assert abs(angle / light.deflection_angle(b, GM_SUN) - 1) < 2e-3
    assert out["max_eps"] < 1e-6 and not out["blocked"]


def test_clock_rates_and_redshift():
    s = two_body()
    r = clocks.body_clock_rates(s.gm, s.pos, s.vel, s.radius)
    assert r[0] < r[1] < 0          # Earth-surface clock is slower than a Moon-surface clock
    assert abs(r[0] + 6.96e-10) < 2e-11
    assert clocks.redshift(r[0], 0.0) > 0      # light climbing out of Earth's well is reddened


def test_fields_tides_hill_roche():
    s = two_body()
    T = fields.tidal_tensor([1e8, 2e7, -3e7], s.gm, s.pos)
    assert abs(np.trace(T)) < 1e-22                    # Laplace: trace-free in vacuum
    w, V = fields.tidal_principal_axes([1e8, 0, 0], [GM_EARTH], [[0, 0, 0]])
    assert abs(w[0] - 2 * GM_EARTH / 1e24) < 1e-18 and abs(abs(V[0, 0]) - 1) < 1e-12
    g = fields.field([[1e8, 0, 0]], [GM_EARTH], [[0, 0, 0]])
    assert abs(g[0, 0] + GM_EARTH / 1e16) < 1e-12
    rep = fields.pair_report(s.bodies, (0, 1))
    assert abs(rep["hill_radius"] / 6.15e7 - 1) < 0.01            # Moon's Hill radius ~61,500 km
    assert abs(rep["roche"]["fluid"] / 1.84e7 - 1) < 0.01         # ~18,400 km
    assert abs(rep["rs_big"] - 8.87e-3) < 1e-5                    # Earth's Schwarzschild radius 8.87 mm
    # the reported L1 is an equilibrium of the co-rotating effective potential
    L1 = np.array(rep["lagrange"]["L1"])
    h = 1e3
    grad = [(fields.effective_potential_pair(L1 + d, s.gm, s.pos, s.vel, (0, 1))
             - fields.effective_potential_pair(L1 - d, s.gm, s.pos, s.vel, (0, 1))) / (2 * h)
            for d in np.eye(3) * h]
    assert np.linalg.norm(grad) < 1e-9 * GM_EARTH / 3.8e8 ** 2 * 1e3
    lines = fields.field_lines(s.gm, s.pos, [[2e8, 1e8, 0]], 2e6, 400, s.radius)
    assert len(lines[0]) > 10
    sheet = fields.potential_sheet(s.gm, s.pos, [0, 0, 0], 5e8, n=32)
    assert sheet["phi"].shape == (32, 32) and sheet["phi"].min() >= sheet["floor"]
