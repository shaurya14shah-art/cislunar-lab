"""Equations of motion of the Circular Restricted Three-Body Problem (CR3BP).

Setting
-------
Two massive bodies ("primaries", e.g. Earth and Moon) move on circles about their
common centre of mass.  A third body (the spacecraft) is so light that it does
not disturb them.  We describe the spacecraft in a frame that *rotates with the
primaries*, so both primaries sit still on the x axis.

Units (non-dimensional)
-----------------------
* length unit  = distance between the primaries
* mass unit    = total mass of the primaries
* time unit    = 1 / (angular speed of the primaries)  -> one revolution = 2*pi
With these units G = 1 and the only parameter left is

    mu = m2 / (m1 + m2)        (Earth-Moon: about 0.01215)

The big primary sits at (-mu, 0, 0) and the small one at (1-mu, 0, 0).

Equations (Szebehely 1967; Koon, Lo, Marsden & Ross, Eq. 2.3.8 with Omega = -U_bar)

    x'' - 2 y' = dOmega/dx
    y'' + 2 x' = dOmega/dy
    z''        = dOmega/dz

    Omega(x,y,z) = (x^2 + y^2)/2 + (1-mu)/r1 + mu/r2 + mu(1-mu)/2

* (x^2+y^2)/2 is the centrifugal potential (the frame spins);
* (1-mu)/r1 and mu/r2 are the two gravitational potentials;
* the constant mu(1-mu)/2 changes no force - it is the customary convention that
  makes the Jacobi constant equal exactly 3 at L4/L5;
* the terms -2y' and +2x' are the Coriolis acceleration.

Jacobi constant:  C = 2*Omega - (x'^2 + y'^2 + z'^2).  It is the one conserved
quantity of the CR3BP and we use its drift as the accuracy monitor everywhere.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp

# Tolerances requested for all "paper-grade" CR3BP propagation.
RTOL = 3e-14
ATOL = 1e-16


def _r1r2(x, y, z, mu):
    r1 = np.sqrt((x + mu) ** 2 + y * y + z * z)       # distance to the big primary
    r2 = np.sqrt((x - 1 + mu) ** 2 + y * y + z * z)   # distance to the small primary
    return r1, r2


def omega(x, y, z, mu):
    """Effective potential Omega (works on scalars or NumPy arrays)."""
    r1, r2 = _r1r2(x, y, z, mu)
    return 0.5 * (x * x + y * y) + (1 - mu) / r1 + mu / r2 + 0.5 * mu * (1 - mu)


def grad_omega(x, y, z, mu):
    """Gradient of Omega = the 'force per unit mass' felt in the rotating frame
    by a body at rest there (gravity + centrifugal)."""
    r1, r2 = _r1r2(x, y, z, mu)
    a, b = (1 - mu) / r1**3, mu / r2**3
    return np.array([x - a * (x + mu) - b * (x - 1 + mu),
                     y - a * y - b * y,
                     -a * z - b * z])


def hessian_omega(x, y, z, mu):
    """3x3 matrix of second derivatives of Omega.

    Physically this is (minus) the tidal tensor plus the centrifugal term: it
    says how the force changes when you step a little away from (x,y,z).  It is
    what makes neighbouring trajectories separate, so it drives both the
    stability analysis and the state transition matrix.
    """
    r1, r2 = _r1r2(x, y, z, mu)
    a3, b3 = (1 - mu) / r1**3, mu / r2**3
    a5, b5 = 3 * (1 - mu) / r1**5, 3 * mu / r2**5
    X1, X2 = x + mu, x - 1 + mu
    H = np.empty((3, 3))
    H[0, 0] = 1 - a3 - b3 + a5 * X1 * X1 + b5 * X2 * X2
    H[1, 1] = 1 - a3 - b3 + (a5 + b5) * y * y
    H[2, 2] = -a3 - b3 + (a5 + b5) * z * z
    H[0, 1] = H[1, 0] = (a5 * X1 + b5 * X2) * y
    H[0, 2] = H[2, 0] = (a5 * X1 + b5 * X2) * z
    H[1, 2] = H[2, 1] = (a5 + b5) * y * z
    return H


def jacobi(state, mu):
    """Jacobi constant C = 2*Omega - v^2.  ``state`` may be (6,) or (6, N)."""
    s = np.asarray(state, dtype=float)
    return 2 * omega(s[0], s[1], s[2], mu) - (s[3] ** 2 + s[4] ** 2 + s[5] ** 2)


def jacobi_jpl(state, mu):
    """Jacobi constant in NASA/JPL's convention, which omits the constant
    mu(1-mu) from Omega:  C_JPL = C - mu(1-mu).  (Found the hard way: our halos
    matched JPL's catalogue in every state component and in period while C
    differed by exactly 0.0120029... = mu(1-mu).)"""
    return jacobi(state, mu) - mu * (1 - mu)


def rhs(t, s, mu):
    """Time derivative of the 6-state [x, y, z, vx, vy, vz]."""
    x, y, z, vx, vy, vz = s
    gx, gy, gz = grad_omega(x, y, z, mu)
    return np.array([vx, vy, vz, 2 * vy + gx, -2 * vx + gy, gz])


def jacobian(s, mu):
    """6x6 matrix A = d(rhs)/d(state).

        A = [[ 0 , I ],
             [ H , K ]]      H = Hessian of Omega,  K = Coriolis block [[0,2,0],[-2,0,0],[0,0,0]]

    Its eigenvalues at an equilibrium point decide whether small displacements
    grow (real positive eigenvalue), oscillate (imaginary pair) or decay.
    """
    A = np.zeros((6, 6))
    A[:3, 3:] = np.eye(3)
    A[3:, :3] = hessian_omega(s[0], s[1], s[2], mu)
    A[3, 4] = 2.0
    A[4, 3] = -2.0
    return A


def rhs_stm(t, y, mu):
    """42-state right-hand side: the 6 state components plus the 36 entries of
    the State Transition Matrix (STM) Phi(t), which obeys  dPhi/dt = A(t) Phi.

    Phi(t) answers: "if I nudge the *initial* state by a tiny d0, how is the
    state at time t changed?"  ->  d(t) = Phi(t) d0.  That is the tool behind
    differential correction, stability indices and invariant manifolds.
    """
    s = y[:6]
    Phi = y[6:].reshape(6, 6)
    out = np.empty(42)
    out[:6] = rhs(t, s, mu)
    out[6:] = (jacobian(s, mu) @ Phi).ravel()
    return out


def propagate(state0, t_span, mu, events=None, dense=False, t_eval=None,
              rtol=RTOL, atol=ATOL, max_step=np.inf):
    """Integrate the CR3BP with the 8th-order Dormand-Prince method (DOP853).

    ``t_span`` may run backwards (t1 < t0); that is how stable manifolds are grown.
    Returns SciPy's solution object (``.t``, ``.y``, ``.t_events`` ...).
    """
    return solve_ivp(rhs, t_span, np.asarray(state0, float), method="DOP853", args=(mu,),
                     rtol=rtol, atol=atol, events=events, dense_output=dense, t_eval=t_eval,
                     max_step=max_step)


def propagate_stm(state0, t_span, mu, events=None, rtol=RTOL, atol=ATOL, dense=False):
    """Integrate state + STM together.  Returns (solution, Phi_final)."""
    y0 = np.concatenate([np.asarray(state0, float), np.eye(6).ravel()])
    sol = solve_ivp(rhs_stm, t_span, y0, method="DOP853", args=(mu,), rtol=rtol, atol=atol,
                    events=events, dense_output=dense)
    return sol, sol.y[6:, -1].reshape(6, 6)


def jacobi_drift(sol, mu):
    """Largest relative change of the Jacobi constant along a solution."""
    C = jacobi(sol.y[:6], mu)
    return float(np.max(np.abs(C - C[0])) / abs(C[0]))
