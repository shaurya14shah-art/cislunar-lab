"""Unit tests for the CR3BP building blocks."""
import numpy as np
import pytest
from cislab.cr3bp import (rhs, jacobian, propagate, propagate_stm, lagrange_points, classify_point,
                          grad_omega, omega, hessian_omega)
from cislab.cr3bp import fastprop as fp
from cislab.cr3bp.lagrange import MU_ROUTH
from cislab.cr3bp.hill import critical_jacobi, hill_slice, regime
from cislab.cr3bp.richardson import Richardson
from cislab.cr3bp.periodic import halo_from_richardson, lyapunov_orbit
from cislab.cr3bp.families import member_with_z
from cislab.cr3bp.manifolds import (eigen_directions, manifold, branch_towards_secondary,
                                    convergence_check, section_curve)
from cislab.cr3bp.poincare import crossings

S = np.array([1.05, 0.02, 0.03, 0.01, 0.2, -0.05])


def test_gradient_and_hessian_match_finite_differences(mu):
    h = 1e-6
    g = grad_omega(*S[:3], mu)
    H = hessian_omega(*S[:3], mu)
    for i in range(3):
        d = np.zeros(3); d[i] = h
        num = (omega(*(S[:3] + d), mu) - omega(*(S[:3] - d), mu)) / (2 * h)
        assert abs(num - g[i]) < 1e-8
        numH = (grad_omega(*(S[:3] + d), mu) - grad_omega(*(S[:3] - d), mu)) / (2 * h)
        assert np.allclose(numH, H[:, i], atol=1e-7)


def test_stm_matches_finite_differences(mu):
    """Phi(t) must equal the numerically differentiated flow."""
    T = 1.5
    _, Phi = propagate_stm(S, (0, T), mu)
    h = 1e-7
    num = np.empty((6, 6))
    for j in range(6):
        d = np.zeros(6); d[j] = h
        num[:, j] = (propagate(S + d, (0, T), mu).y[:, -1] - propagate(S - d, (0, T), mu).y[:, -1]) / (2 * h)
    assert np.allclose(Phi, num, rtol=1e-5, atol=1e-5)
    assert abs(np.linalg.det(Phi) - 1) < 1e-8        # Liouville: phase-space volume is conserved


def test_compiled_integrator_matches_scipy(mu):
    ref = propagate(S, (0, 3.0), mu).y[:, -1]
    assert np.max(np.abs(fp.integrate(S, 0, 3.0, mu) - ref)) < 1e-11
    _, P = propagate_stm(S, (0, 3.0), mu)
    _, Pf = fp.integrate_stm(S, 0, 3.0, mu)
    assert np.max(np.abs(P - Pf)) / np.max(np.abs(P)) < 1e-9
    back = fp.integrate(fp.integrate(S, 0, 3.0, mu), 3.0, 0.0, mu)     # there and back again
    assert np.max(np.abs(back - S)) < 1e-10


@pytest.mark.parametrize("m", [1e-6, 3.0e-6, 0.001, 0.01215, 0.03, 0.2, 0.5])
def test_lagrange_points_are_equilibria_for_any_mu(m):
    for name, L in lagrange_points(m).items():
        assert np.linalg.norm(grad_omega(*L, m)) < 1e-11, name


def test_stability_classification(mu):
    for p in ("L1", "L2", "L3"):
        assert classify_point(mu, p)["unstable"]
    for p in ("L4", "L5"):
        assert not classify_point(mu, p)["unstable"]
    assert classify_point(0.05, "L4")["unstable"]            # above the Routh limit
    assert abs(MU_ROUTH - 0.0385208965) < 1e-9
    assert abs(classify_point(mu, "L4")["jacobi"] - 3.0) < 1e-12


def test_hill_region_gateways(mu):
    c = critical_jacobi(mu)
    assert c["L1"] > c["L2"] > c["L3"] > c["L4"]
    # just above C(L1) the neck at L1 is closed (forbidden), just below it is open
    L1 = lagrange_points(mu)["L1"]
    assert 2 * omega(*L1, mu) - (c["L1"] + 1e-6) < 0
    assert 2 * omega(*L1, mu) - (c["L1"] - 1e-6) > 0
    assert "L1 neck open" in regime(mu, 0.5 * (c["L1"] + c["L2"]))
    u, v, W = hill_slice(mu, 3.17, n=50)
    assert W.shape == (50, 50)


def test_richardson_seed_is_close_and_converges(mu):
    R = Richardson(mu, "L2")
    seed = R.state(0.03)
    orb, _ = halo_from_richardson(mu, 0.03)
    assert orb.residual < 1e-11
    assert np.max(np.abs(seed - orb.state0)) < 5e-3           # third-order theory: good to ~1e-3
    assert abs(R.period(0.03 / R.gamma) - orb.period) / orb.period < 5e-3
    # the orbit really is periodic: fly a full period and come back
    end = fp.integrate(orb.state0, 0, orb.period, mu)
    assert np.max(np.abs(end - orb.state0)) < 1e-9
    o1, _ = halo_from_richardson(mu, 0.02, "L1")
    assert o1.residual < 1e-11 and o1.state0[0] < 1 - mu


def test_lyapunov_orbit_is_planar_and_periodic(mu):
    orb, _ = lyapunov_orbit(mu, 0.02)
    end = fp.integrate(orb.state0, 0, orb.period, mu)
    assert np.max(np.abs(end - orb.state0)) < 1e-9 and orb.state0[2] == 0


def test_monodromy_structure(mu, family):
    """Multipliers come in reciprocal pairs, two equal 1, determinant 1."""
    orb = member_with_z(family, 0.08)
    lam = np.sort(np.abs(orb.multipliers))
    assert abs(lam[0] * lam[-1] - 1) < 1e-6
    assert np.sum(np.abs(orb.multipliers - 1) < 1e-4) == 2
    assert abs(np.linalg.det(orb.monodromy) - 1) < 1e-6


def test_family_is_monotone_in_stability(family):
    """Along the classic part of the family the halos get steadily less unstable."""
    k = np.array([o.stability_index for o in family])
    z = np.array([abs(o.state0[2]) for o in family])
    first = slice(0, int(np.argmax(z)))
    assert np.all(np.diff(k[first]) < 0)
    assert k[0] > 500 and k.min() < 1.0001


def test_manifolds(mu, family):
    orb = member_with_z(family, 0.08)
    lu, vu, ls, vs = eigen_directions(orb)
    assert abs(lu * ls - 1) < 1e-6 and lu > 1
    # the eigenvector really is one: M v = lambda v
    assert np.linalg.norm(orb.monodromy @ vs - ls * vs) / np.linalg.norm(vs) < 1e-6
    sg = branch_towards_secondary(orb, "stable")
    man = manifold(orb, "stable", 24, 1e-6, sg, section_x=1 - mu)
    ph, cur = section_curve(man)
    assert len(cur) >= 15 and np.allclose(cur[:, 0], 1 - mu, atol=1e-12)
    # stable-manifold property: flown FORWARD from the section, the trajectory
    # arrives back at the halo (distance ~ eps)
    k = len(man) // 3
    while man[k].crossing is None:
        k += 1
    end = fp.integrate(man[k].crossing, man[k].t_cross, 0.0, mu)
    on_orbit = fp.integrate(orb.state0, 0.0, man[k].phase * orb.period, mu)
    assert 1e-9 < np.linalg.norm(end[:3] - on_orbit[:3]) < 2e-6      # <= eps, as constructed
    # Jacobi constant on the manifold equals that of the orbit
    from cislab.cr3bp import jacobi
    assert abs(jacobi(man[k].crossing, mu) - orb.jacobi) < 1e-5
    chk = convergence_check(orb, "stable", 1e-6, sg)
    assert chk["ok"] and chk["relative"] < 0.01, chk


def test_poincare_crossings(mu, family):
    orb = member_with_z(family, 0.08)
    t, s = crossings(orb.state0, mu, 2.2 * orb.period, axis=1, value=0.0, direction=0)
    assert len(t) >= 4 and np.max(np.abs(s[:, 1])) < 1e-12
    assert abs(t[1] - orb.period) < 1e-8       # second y=0 crossing after t=0 is one period later
