# Decisions log

Newest last. Each entry: what was decided, why, and what it costs.

1. **Stack: Python core + FastAPI/WebSocket server + Three.js in the browser.**
   Chosen with the user. Python remains the single source of truth; the browser
   never computes physics. Cost: one command to start a local server.

2. **mu and units from JPL's periodic-orbit API** (mu = 0.01215058560962404,
   L = 389,703.2648 km, T = 382,981.2891 s) rather than the brief's rounded
   0.01215058. Why: it lets published halo orbits be reproduced digit for digit.
   The two differ in the 8th decimal.

3. **Jacobi constant with the mu(1-mu) term** (Szebehely / Koon et al.), so
   C(L4) = 3. Discovered during verification that JPL omits it; `jacobi_jpl()`
   converts. Cost: our C differs from JPL's by 0.012003.

4. **A compiled copy of DOP853 (`fastprop.py`).** The brief asks for DOP853 via
   `solve_ivp`; that path exists (`core.propagate`) and is the reference. The
   compiled copy uses SciPy's own coefficients and is tested equal to it. Why:
   the family continuation was taking minutes and the Monte Carlo would have
   taken days. Cost: two code paths to keep in agreement (one test does that).

5. **Pseudo-arclength continuation** instead of stepping in z0. The L2 family's
   far-side z is not monotonic (it peaks near 79,000 km), so natural-parameter
   continuation would stall there.

6. **"Northern" follows JPL** (z > 0 at the far-side crossing). Richardson's
   class label refers to the near-side crossing and is the opposite for L2.

7. **Injection point = where the Moon-bound stable manifold first crosses the
   Moon plane x = 1 - mu.** A definite, reproducible place that is also the
   natural Poincare section. Cost: results are specific to this choice; an
   injection further from the Moon would see different error growth.

8. **Capture definitions.** A: within d_tol = 5,000 km of the halo for one
   consecutive revolution. B: within two Hill radii of the Moon, no impact,
   until (nominal arrival + one period). Chosen with the user to report both.
   First draft used a 2-period horizon for B and 1 for A, which made the
   "looser" definition stricter in time; fixed so both are judged over the same
   duration. Sensitivity to d_tol and to the number of revolutions is reported.

9. **Unusable phases are excluded and counted.** Manifold trajectories that hit
   the Moon, or whose unperturbed forward flight fails to hold the halo for two
   periods (because they graze the Moon and amplify even round-off), are not
   used as injection points. The count is stored (`scan`) and reported. Cost: a
   selection effect: pooled tolerances describe *usable* routes only.

10. **Monte Carlo tolerances 1e-12 / 1e-14, step cap 0.01**, looser than the
    3e-14 used for orbit construction, for speed. Justified by the halving
    test: outcomes are unchanged.

11. **1PN in the CR3BP: test-particle EIH term with the primaries on their
    Newtonian circle.** The primaries' own 1PN correction only changes the
    relation between separation and period at the 1e-11 level, absorbed by
    using observed units. Nominal trajectories and the reference halo are *not*
    re-targeted in the 1PN model, which is the pessimistic choice.

12. **1PN runs use the adaptive integrator.** EIH forces depend on velocity, so
    the kick-drift splitting behind Yoshida's method no longer applies. The
    server refuses the combination with an explanation.

13. **Collisions flagged, not merged; spin displayed, not dynamical.** Both
    would need physics this code does not contain.

14. **Regime bands are conventions** (boundaries at 1e-8, 1e-3, 3e-2, 1e-1 of
    max(eps, beta^2)), stated as such; raw numbers are always shown.

15. **Interface direction: "navigation chart".** Two directions were
    considered. (a) *Chart*: pale water for allowed space, buff land for
    forbidden regions, ink lines, slanted serif for feature names as on sea
    charts, magenta reserved for overlays, here everything relativistic.
    (b) *Deep-field instrument*: dark scene with luminous trails and glass
    labels. Chosen (a): it turns the hardest idea (energy decides where you can
    go) into something read at a glance, it prints well for the paper, and it
    is not the dark-dashboard look the brief rules out.

16. **Controls are contextual.** Notes open beside the object clicked; the left
    panel is a chart key whose rows show the real line styles and toggle
    layers; insets hold secondary plots. No tab bar, no card grid.

17. **Zero-velocity surface drawn as stacked contour lines**, not a solid
    mesh: it reads as a chart, keeps the interior visible, and needs no
    marching-cubes dependency.

18. **Light rays in the N-body scene use the weak-field model** and are drawn
    dashed where it is not trustworthy; exact rays are in the single-body
    close-up. One honest approximation plus one exact special case, rather than
    a multi-body strong-field ray tracer that would be beyond this project.

19. **"Fields" means gravitational fields.** Magnetic dipole fields were
    considered as a stretch goal and not implemented: they do not act on a
    neutral spacecraft's trajectory, so they would be decoration here.

20. **Export: PNG and CSV.** (The brief's "STEP-free" is read as "no CAD export
    needed".)

21. **Manifold starts use constant Floquet amplitude, not constant step.**
    *This entry records a mistake and its correction.* The first full Monte
    Carlo used the textbook recipe: at evenly spaced phases, step a fixed eps
    along the local stable eigenvector. The robustness run with eps/10 then
    moved the pooled 50% position tolerance far outside its bootstrap range
    (198 m -> 79 m at z0 = 0.04; 530 m -> 284 m at 0.10; 664 m -> 2,480 m at
    0.16), although the tube had not moved (section-curve convergence 3e-5).
    Diagnosis: a trajectory on the tube is identified by its Floquet label
    theta = phase - ln(eps/|p(phase)|)/ln(lambda), where p is the Floquet-periodic
    eigenvector. Constant eps makes theta a non-uniform function of phase, so
    trajectories are sampled unevenly, and changing eps slides that uneven
    pattern around. Tolerance varies by a factor of several hundred between
    trajectories, so the pooled value changed. The plateau-then-jump in the
    first run's tolerance-versus-amplitude curve was the same artefact.
    Fix: step eps*|p(phase)|/max|p| (`manifold(..., uniform=True)`, now the
    default). Then theta = phase + constant: even sampling, and eps only slides
    it. Check: linear-theory pooled tolerances over all usable routes agree
    between eps and eps/10 to 0.1% for large halos and ~10% for the smallest.
    The whole experiment was rerun. The first run's summary is kept as
    `data/mc_summary_v1_uneven_sampling.json` for the record.

22. **Monte Carlo budget moved from draws to routes: 24 routes x 150 draws**
    (was 12 x 300) per error level. Route-to-route variation dominates the
    uncertainty, so this is the better use of the same compute.

23. **Bootstrap ranges resample routes as well as draws**, but with only 24
    routes they still under-represent the heavy tail of route sensitivities.
    `scripts/route_study.py` therefore also evaluates the (Monte-Carlo-
    validated) linear prediction on every usable route of a 96-trajectory scan.
