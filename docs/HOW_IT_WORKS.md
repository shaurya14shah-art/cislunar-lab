# How it works

> New to this? Read `docs/PLAIN_GUIDE.md` first (everyday language, every graph explained) and keep `docs/GLOSSARY.md` open. This file is the fuller, more technical version.

A plain-language guide to the physics and the code, written so you can explain
every part to a teacher or a judge. Each section says what the idea is, why it
is true, where it lives in the code, and which test proves the code does it.

Formulas are written in plain text. `x'` means "rate of change of x", `x''`
means acceleration.

---

## 1. The stage: three bodies, one of them negligible

**Idea.** Earth and Moon orbit their common centre of mass on (nearly) circles.
A spacecraft is so light that it does not disturb them. That is the *Circular
Restricted Three-Body Problem* (CR3BP).

**Trick 1: stand on a turntable.** Watch from a frame that rotates once a month
with the Earth-Moon line. Now Earth and Moon stand still, and the only thing
moving is the spacecraft. The price: two "fictitious" forces appear, exactly as
on a merry-go-round.

* centrifugal force, pushing outward, strength = distance from the axis;
* Coriolis force, pushing sideways on anything that moves, strength = 2 x speed.

**Trick 2: choose units that make the numbers disappear.**

| Quantity | Unit | Earth-Moon value |
|---|---|---|
| length | Earth-Moon distance | 389,703 km |
| time | 1 / (angular speed) | 382,981 s = 4.43 days (so one month = 2 pi) |
| mass | Earth + Moon | |

With these units G = 1 and a single number is left: `mu = m_Moon / (m_Earth + m_Moon) = 0.01215`.
Earth sits at x = -mu, the Moon at x = 1 - mu.

**Equations of motion** (`cislab/cr3bp/core.py`, function `rhs`):

    x'' - 2 y' = dOmega/dx
    y'' + 2 x' = dOmega/dy
    z''        = dOmega/dz

    Omega = (x^2 + y^2)/2  +  (1 - mu)/r1  +  mu/r2  +  mu(1 - mu)/2

* `(x^2+y^2)/2`: centrifugal potential.
* `(1-mu)/r1`, `mu/r2`: gravity of Earth and Moon (r1, r2 are the distances to them).
* `-2y'`, `+2x'`: Coriolis.
* the last constant changes no force; it is a convention (see section 3).

**Where it is derived.** Koon, Lo, Marsden & Ross, *Dynamical Systems, the
Three-Body Problem and Space Mission Design*, Eq. (2.3.8); Szebehely, *Theory of
Orbits* (1967). To derive it yourself: write Newton's law in the inertial frame,
substitute `X = x cos t - y sin t`, `Y = x sin t + y cos t`, and collect terms.

**Tests.** `test_01_mu_zero_is_kepler` switches the Moon off (mu = 0) and checks
that a circular orbit stays circular and an ellipse closes after exactly one
Kepler period. If the Coriolis or centrifugal signs were wrong, this would fail.

---

## 2. Lagrange points: where you can park

**Idea.** A Lagrange point is where a spacecraft at rest (in the rotating frame)
stays at rest: gravity of both bodies and the centrifugal force cancel. In
symbols, the slope of Omega is zero.

* **L1, L2, L3** lie on the Earth-Moon line. Their positions solve a
  fifth-degree equation with no formula, so the code finds the root numerically
  (`cislab/cr3bp/lagrange.py`, Brent's method).
* **L4, L5** make equilateral triangles with Earth and Moon. That is exact.

**Stable or not?** Nudge the spacecraft slightly and ask how the nudge evolves.
Near the point the equations become linear: `d' = A d`, where A is a 6x6 matrix
(the *Jacobian*) built from the second derivatives of Omega plus the Coriolis
terms. The solutions are `e^(s t)` with s an eigenvalue of A.

| Eigenvalue s | Meaning |
|---|---|
| real, positive | grows exponentially: unstable |
| real, negative | decays |
| purely imaginary (+- i w) | oscillates forever with frequency w |

At L1 and L2 you get one growing, one decaying, and two oscillating pairs
("saddle x centre x centre"). At L2: growth rate lambda = 2.1587, so an error
multiplies by e every 1/2.1587 time units = 2.05 days. The two oscillations
(frequencies nu = 1.8626 in the plane, omega_z = 1.7862 out of it) are the raw
material of halo orbits.

At L4/L5 all eigenvalues are imaginary provided mu < 0.03852 (Routh's limit), so
Earth-Moon L4/L5 are stable.

**Tests.** `test_03` (positions to 6 figures, and to 13 against JPL), `test_04`
(the three rates, *and* an experiment: displace a particle by 1e-8 at L1,
integrate the full nonlinear equations, fit the growth rate, compare with
lambda), `test_02` (a particle at L4 stays for 100 time units).

---

## 3. The Jacobi constant and forbidden regions

**Idea.** Energy is not conserved in a rotating frame, but one combination is:

    C = 2 Omega - (speed)^2          "Jacobi constant"

High C = low energy. Since speed^2 cannot be negative, a spacecraft with a given
C can only be where `2 Omega >= C`. The boundary is the *zero-velocity surface*;
the allowed side is the *Hill region*.

Lower C a little at a time and gates open in order: first at L1 (Earth-Moon
traffic becomes possible), then L2 (escape past the Moon), then L3, and at C = 3
nothing in the plane is forbidden any more. In the interface this is the "sea
level" slider: land is forbidden, water is allowed.

**Why we care numerically.** C must stay constant along any computed
trajectory. Its drift is therefore a free accuracy meter, shown at all times in
the interface. `test_05` requires drift < 1e-11.

**A trap we fell into.** NASA/JPL's catalogue defines Omega *without* the
constant mu(1-mu)/2, so its C is smaller than ours by mu(1-mu) = 0.012003. Our
halos matched JPL in every coordinate and in period while C disagreed by exactly
that amount. Both conventions are legitimate; `jacobi_jpl()` converts.

Code: `cislab/cr3bp/hill.py`.

---

## 4. Halo orbits: closing the loop

**Idea.** Near L2 the in-plane and out-of-plane oscillations have slightly
different frequencies (1.8626 vs 1.7862), so a small orbit never closes. Make
the orbit *bigger* and nonlinear effects shift the frequencies; at one
particular size they become equal and the path closes into a tilted loop: a
*halo orbit*. From the Earth it looks like a halo around the Moon.

### 4a. A first guess: Richardson's formula

Richardson (1980) expanded the equations as a power series in the orbit size
(to third order) and adjusted the frequency as he went so the series does not
blow up (the Lindstedt-Poincare method). The result is an explicit formula
`x(t), y(t), z(t)` made of cosines and sines. It is only approximate: for a
7,800 km halo it is off by a few hundred kilometres. Good enough to *start* from.
Code: `cislab/cr3bp/richardson.py`.

### 4b. Making it exact: differential correction

**Symmetry trick.** If a path leaves the x-z plane at right angles and later
hits that plane at right angles again, its mirror image completes a closed
orbit. So we need: at the second crossing, `vx = 0` and `vz = 0`.

**Shooting.** Guess the start `(x0, z0, vy0)` and the half-period `tau`. Fly.
The crossing will not be perpendicular; we get errors. Newton's method asks:
*how would the end point change if I changed the start?* That is exactly what
the next tool provides.

### 4c. The state transition matrix (STM)

The STM `Phi(t)` is a 6x6 matrix that answers: "if I change the initial state by
a tiny d0, the state at time t changes by Phi(t) d0." It obeys

    Phi' = A(t) Phi,     Phi(0) = identity,

where A(t) is the same Jacobian as in section 2, evaluated along the path. So we
integrate 6 + 36 = **42 numbers** together (`rhs_stm`). With Phi in hand, each
Newton step is one small linear solve. Typically 3-5 steps reach 1e-13.

`test_stm_matches_finite_differences` checks Phi against brute-force
differencing.

### 4d. Walking along the family: continuation

Three conditions, four unknowns, so halos come in a one-parameter *family*. To
step to the next member we predict along the family's tangent and correct with
the extra condition "stay on the plane perpendicular to the tangent"
(*pseudo-arclength continuation*). Unlike "increase z0 a bit each time", this
keeps working when the family curve turns around, which it does: the far-side z
peaks at about 79,000 km and then decreases towards the Moon.
Code: `cislab/cr3bp/periodic.py`, `families.py`.

### 4e. How unstable is a given halo?

Integrate the STM once around the orbit: that is the *monodromy matrix* M. Its
eigenvalues ("multipliers") say what one revolution does to an error. They come
in reciprocal pairs; the largest, lambda_u, is the factor by which an error
grows per revolution. The **stability index** is `k = (lambda_u + 1/lambda_u)/2`.

| Halo (far-side z) | k | error growth per revolution |
|---|---|---|
| 7,800 km | 588 | x 1176 |
| 39,000 km | 283 | x 565 |
| 70,000 km | 32 | x 64 |
| near-rectilinear (NRHO) | about 1 | none |

**Tests.** `test_06_published_halo` rebuilds three orbits from JPL's catalogue
starting from *our own* Richardson seed and continuation, taking only z0 from
JPL; x0, vy0, period and C agree to better than 1 part in 10 million and k to
1 part in 10,000.

---

## 5. Invariant manifolds: free rides

**Idea.** An unstable orbit is like a pencil balanced on its tip, but with a
twist: there is one special direction along which things fall *towards* it.

* Step off along the **unstable eigenvector** of M and you drift away: all such
  paths form the *unstable manifold*.
* Step off along the **stable eigenvector** and integrate *backwards in time*:
  you trace the paths that, flown forwards, arrive on the orbit for free. That
  is the *stable manifold*, the basis of ballistic capture.

Do this from every point of the orbit (transporting the eigenvector with the
STM) and you get a tube in space.

**The step epsilon.** We start at `orbit + eps * v` with eps = 1e-6 length
units (about 390 m). Smaller eps = more accurate but the path lingers longer.
`convergence_check` repeats with eps/10 and measures how far the tube moved
where it crosses the Moon plane.

**A sampling pitfall we hit.** "Step eps off the orbit at evenly spaced
phases" draws the correct tube but does *not* pick its trajectories evenly: how
far along the tube a start lands depends on how large the eigenvector is at
that phase. Worse, changing eps reshuffles which trajectories get picked. Our
first experiment's pooled answer changed by up to a factor of 4 when eps
changed, and that is how we found it. The cure is to scale each step by the
size of the *Floquet-periodic* eigenvector there (the eigenvector with its
overall exponential growth divided out). Every start then has the same
"Floquet amplitude", the starts are evenly spread over the trajectories, and
eps only slides them along. This is `uniform=True`, the default.

**Poincare section and capture curve.** Cut the tube with the plane through the
Moon, `x = 1 - mu`. The cut is a closed loop: the *capture curve*. A spacecraft
crossing that plane on the loop is delivered to the halo.
Code: `cislab/cr3bp/manifolds.py`, `poincare.py`.

**Tests.** `test_manifolds`: the stable eigenvector really satisfies M v = lambda v;
a stable-manifold state flown forward returns to within eps of the orbit; the
Jacobi constant on the tube equals the orbit's; the eps and eps/10 curves agree.

---

## 6. The experiment: how big an error can capture absorb?

This is the paper's engine (`cislab/cr3bp/montecarlo.py`).

1. **Nominal.** For each halo, take 24 points on the capture curve (spread
   evenly around the tube, skipping the ones that hit or graze the Moon).
2. **Perturb.** Add a random error: position error drawn from a bell curve of
   width sigma_r per axis, or velocity error of width sigma_v.
3. **Coast.** No correction burns. Integrate.
4. **Classify.**
   * **A, halo-tracking:** stays within 5,000 km of the halo for one full revolution.
   * **B, Moon-region:** stays within two Hill radii (124,000 km) of the Moon,
     without impact, for the same length of time.
   * otherwise *impacted* (hit Moon or Earth) or *escaped*.
5. **Repeat** for 150 draws x 24 points x 22 error sizes x 10 halos, and read off
   the error size at which half are captured (the "50% tolerance").

**Why two definitions?** Because the answer depends on what "captured" means,
and a paper that hides that would be overclaiming. B is far more forgiving: a
spacecraft that misses the halo but falls towards the Moon still counts.

**Common random numbers.** The same 150 random draws are reused for every error
size, halo and model. Then when two curves differ, it is because of the thing
that changed, not luck.

**A check from theory.** The STM predicts the answer without any Monte Carlo:
an error is stretched by a factor g (largest singular value of the position
part of the STM from injection to one revolution after arrival). Capture needs
`g * error < 5000 km`, so the 50% tolerance should be about `5000 km / (0.674 g)`.
Route by route, the simulated tolerance is between 0.6 and 4.8 times this
prediction across four orders of magnitude of g (paper, Fig. 7); pooled over
routes the two agree within a factor of about 2.

**Speed.** The experiment needs millions of trajectories. `fastprop.py` is the
same DOP853 integrator as SciPy's (same coefficients, imported from SciPy) but
compiled with Numba; `test_compiled_integrator_matches_scipy` holds them equal
to 1e-11.

---

## 7. Integrators: why two kinds

| | DOP853 (adaptive Runge-Kutta, 8th order) | Yoshida (symplectic, 4th order) |
|---|---|---|
| step size | chosen automatically from an error estimate | fixed |
| energy error | tiny, but slowly accumulates | larger, but oscillates and never drifts |
| good for | close encounters, high accuracy, velocity-dependent forces | very long N-body runs |

A symplectic method takes alternating "drift" (move with current velocity) and
"kick" (change velocity with current force) sub-steps. Yoshida's recipe chains
three such steps with weights (w1, w0, w1); w0 is *negative*, and that backward
middle step is what cancels the leading error. `test_yoshida_is_fourth_order`
halves the step and checks that the error falls 16-fold.

---

## 8. The N-body sandbox

`cislab/nbody/`. Any number of point masses under Newton's law
`a_i = sum_j G M_j (r_j - r_i)/|r_j - r_i|^3`.

* **Frames.** Inertial; barycentric (centre of mass at rest); co-rotating with
  any pair (origin at the pair's barycentre, x axis from one to the other).
* **Collisions** are *flagged*, never merged: merging needs physics this code
  does not have, and silently merging would hide that.
* **Spin** is stored and displayed but does not affect motion (point masses).
* **Units.** SI internally. `unit_system()` prints the conversion to pair units,
  which is exactly how the CR3BP units of section 1 arise.

---

## 9. Relativity, and when to believe it

Two small numbers measure how relativistic a situation is:

    eps  = G M / (r c^2)      how deep in a gravity well
    beta^2 = (v/c)^2          how fast

For orbits the two are similar in size. At the Moon's distance from Earth,
eps = 1.15e-11.

### 9a. Three levels of theory

| Theory | Keeps terms of size | Throws away | Valid when |
|---|---|---|---|
| Newton | 1 | eps | eps tiny |
| 1PN (Einstein-Infeld-Hoffmann) | 1 and eps | eps^2 | eps below about 0.01 |
| Exact (Schwarzschild) | everything | nothing, but only one massive body | any eps |

The **regime meter** shows eps and beta^2 for every body live and names the
band. It is honest about the bands being conventions, and shows the raw numbers.

### 9b. 1PN N-body dynamics

`cislab/nbody/kernels.py`, `accel_eih_correction`, documents each term. In
words, compared with Newton:

* a body deep in other bodies' potential feels a modified pull (the `-4 sum mu/r` term);
* the source's own potential energy and kinetic energy gravitate too;
* there are velocity-dependent forces with no Newtonian analogue;
* a source that is itself accelerating pulls differently.

**Tests.** Mercury's perihelion advance comes out at 42.980 arcseconds per
century (`test_07`), and exactly 0 with the term switched off. The correction
scales exactly as 1/c^2. The 1PN energy is conserved hundreds of times better
than the Newtonian energy under 1PN dynamics.

### 9c. Schwarzschild: exact motion around one mass

`cislab/relativity/schwarzschild.py`. Radial motion obeys

    r'' = - M/r^2 + L^2/r^3 - 3 M L^2 / r^4        (G = c = 1, ' = per proper time)

The first two terms are Newton's gravity and the centrifugal barrier. The third
is new, attractive, and wins at small r. Consequences, all tested (`test_08`):

* **precession**: orbits do not close;
* **ISCO at r = 6 GM/c^2**: closer in, no circular orbit is stable. Found
  numerically as the radius of minimum angular momentum: 6.0000000003;
* **photon sphere at 3 GM/c^2**: light can circle;
* **light bending**: 1.7512 arcseconds at the Sun's edge.

### 9d. Clocks

A clock at depth eps moving at speed v ticks at `1 - eps - beta^2/2` relative to
a far-away clock at rest. For GPS: the satellite is higher (ticks faster by
45.7 microseconds/day) and moving (slower by 7.1), net **+38.58** (`test_09`).
In the interface each body carries a dial; because real differences are parts
per billion, the *lag* is magnified by a stated factor.

### 9e. Does relativity matter for the paper's question?

The 1PN acceleration on a halo spacecraft is about 4e-11 of the Newtonian one.
Three independent measurements in the experiment say it does not change the
capture tolerance (paper, section 4.5): identical capture curves, 11 of about
3.2 million outcome classifications changed, and an "equivalent injection
error" of a few millimetres against tolerances of metres to kilometres.

---

## 10. Fields

`cislab/fields/gravity.py`.

* **Potential** `Phi = - sum GM/r`, drawn as a rubber sheet (deep wells are cut
  off at the deepest 1% so they do not become infinite spikes; the display says so).
* **Field lines** follow the direction of gravity.
* **Tidal tensor** `T = sum GM (3 n n^T - I)/r^3`: tells how two nearby particles
  are pulled apart (+2GM/r^3 along the line to the mass) or squeezed (-GM/r^3
  across it). Its three eigenvalues always add to zero.
* **Hill sphere** `a (m/3M)^(1/3)`: where the small body's gravity beats the big
  body's tide. **Roche limit**: where the big body's tide beats the small body's
  self-gravity.
* **Lagrange points of any pair**: the CR3BP points, scaled and rotated to the
  pair; flagged when the pair's orbit is not circular.

---

## 11. Where each claim is tested

| Claim | Test |
|---|---|
| Equations of motion | `test_01`, `test_gradient_and_hessian_match_finite_differences` |
| Lagrange points and stability | `test_02`, `test_03`, `test_04`, `test_stability_classification` |
| Jacobi conservation | `test_05` |
| STM | `test_stm_matches_finite_differences` |
| Halo orbits | `test_06`, `test_richardson_seed_is_close_and_converges` |
| Manifolds | `test_manifolds` |
| Monte Carlo robustness | `test_10`, plus the batch robustness runs |
| Integrators | `test_yoshida_is_fourth_order`, `test_compiled_integrator_matches_scipy` |
| 1PN | `test_07`, `test_eih_reduces_to_newton_and_scales_as_c_minus_2`, `test_test_particle_1pn_matches_general_eih` |
| Schwarzschild, light, clocks | `test_08`, `test_09`, `test_schwarzschild_orbits_and_clocks` |
| Fields | `test_fields_tides_hill_roche` |
| Server | `tests/test_server.py` |

## 12. Questions a judge might ask

**"Why trust a simulation?"** Because it reproduces ten independent known
results, several to 10+ digits, and because halving every numerical tolerance
changes none of the conclusions.

**"Is this a real mission design?"** No. The CR3BP leaves out the Sun (which
perturbs Earth-Moon halos noticeably over a month), the Moon's orbital
eccentricity (0.055), and the Moon's lumpy gravity. Real missions also fly
correction burns every few days. What the study gives is how *sensitivity*
scales with halo size, which is a property of the dynamics.

**"Why is the tolerance so small?"** Because nothing corrects the error. An
unstable orbit multiplies errors by hundreds per revolution. The result is best
read as "how soon must the first correction burn come, and how precise must
navigation be until then".

**"What did relativity change?"** Nothing measurable, and the code shows why
with numbers rather than assertion.
