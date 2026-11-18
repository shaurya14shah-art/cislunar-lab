# Final summary

## What is done

| Deliverable | Where | State |
|---|---|---|
| Physics core (CR3BP, N-body, 1PN, Schwarzschild, clocks, light, fields) | `cislab/` | complete, tested |
| Local server (REST + WebSocket) | `cislab/server/` | complete, tested |
| Interactive 3D front end, two sheets | `web/` | complete; checked by automated screenshots |
| Monte Carlo experiment, both capture definitions, Newtonian and 1PN | `scripts/run_experiment.py`, `data/mc/` | complete: 3.17 million trajectories in the main run plus robustness and all-routes runs |
| Paper draft with 8 figures and 3 tables | `paper/paper.tex`, `paper/paper.pdf` | complete draft, every number generated from data |
| One command that regenerates everything | `scripts/reproduce_all.py` | complete |
| Verification report | `docs/VERIFICATION.md` | 17 of 17 checks pass; 50 of 50 tests pass |
| How-it-works guide, sources, decisions log | `docs/` | complete |

## Headline results (Earth-Moon CR3BP, uncorrected coast, injection at the Moon plane)

* Holding the halo for one revolution (definition A): the position error at which half the spacecraft are
  captured rises from about **100 m** (smallest halo, Az = 7,800 km, k = 588) to about **18 km** (largest, Az = 76,000 km, k = 11);
  in velocity from **4 mm/s** to **185 mm/s**. Against the stability index: roughly k^-1.22 (position) and k^-0.95 (velocity).
* Merely staying within two Hill radii of the Moon (definition B) tolerates errors 60 to 140 times larger; for the two
  largest halos more than half were still "captured" at the largest error tested, so only a lower bound exists.
* The arrival route matters more than the halo: tolerance differs by factors of 15 to about 1,000 between routes of one halo.
  Routes passing close to the Moon are the fragile ones.
* A linear prediction from the state transition matrix matches the simulation route by route across four orders of magnitude.
* Each extra half revolution of required dwell costs a factor of about sqrt(2k); tolerance scales as the square of the capture-tube width.
* **Relativity:** the 1PN acceleration is 4e-11 of the Newtonian one. With identical random draws, 11 of 3,168,000 outcome
  classifications changed. Relativistic corrections do not change the capture tolerance.

## What is verified

All ten checks of the brief, by automated tests whose measured values are in `docs/VERIFICATION.md`. Highlights:
L1/L2 to 13 digits against JPL; three published JPL halo orbits (including a near-rectilinear one) reproduced to
better than 1 part in 10^11 from the code's own seed; Mercury 42.980 arcsec/century; ISCO 6.000000000 GM/c^2;
solar-limb deflection 1.7512 arcsec; GPS +38.575 microseconds/day; Jacobi drift below 1e-13; halved tolerances change
0.04 to 0.08% of individual outcomes and no conclusion.

## A mistake found and fixed (worth telling a judge)

The first full run used the standard recipe for manifolds (a fixed step eps at evenly spaced phases). The robustness
test then failed: eps/10 moved the pooled tolerance by up to a factor of 4. Cause: that recipe samples the tube's
trajectories unevenly, and the unevenness shifts with eps. Fix: constant Floquet amplitude. The experiment was rerun;
eps and eps/10 now agree within the error bars, and an artificial plateau in the first run's results disappeared.
Details: `docs/DECISIONS.md` entry 21 and section 4.7 of the paper.

## What is NOT verified, or is limited

* **Sources not read in full.** Richardson (1980) and Howell (1984) are paywalled: only abstracts were opened.
  Richardson's coefficients were cross-checked against the open-source Orekit library. Szebehely (1967), Einstein-Infeld-Hoffmann
  (1938) and Landau-Lifshitz were not opened; the EIH equations were taken from Holman et al. (2023), Eq. 2. The IAU nominal
  solar radius is quoted from general knowledge. `docs/SOURCES.md` lists every case and what protects the code against it.
* **The 1PN energy formula** is validated only by a conservation test, not against a primary source.
* **The model.** No Sun, no lunar eccentricity, no lunar gravity harmonics, no correction burns, isotropic Gaussian errors.
  Absolute tolerances would change in a realistic model; whether the scaling survives is untested.
* **Route selection.** 0 to 29% of scanned manifold trajectories per halo were excluded because they hit or graze the Moon.
  Results describe usable routes.
* **Pooled statistics rest on 24 routes per halo.** The all-routes check (35 to 48 routes) agrees within error bars, but the
  bootstrap ranges still under-represent route-to-route spread.
* **Definition-dependence.** Tolerances are relative to the 5,000 km tube and one-revolution dwell; both choices change the numbers a lot.
* **Front end.** Exercised through automated screenshots of the main interactions on Chromium with software rendering.
  Not tested on Safari or Firefox, on a phone, or on your Mac's GPU; no automated interaction suite.
* **Halo, manifold and Monte Carlo tools in the interface** are set up for the Earth-Moon mass ratio. Lagrange points,
  stability and forbidden regions work for any mass ratio.
* **Spin** is displayed but not dynamical; **collisions** are flagged, not modelled; **light rays** in the N-body scene are
  weak-field (flagged where invalid); exact strong-field motion exists only for a single non-rotating mass.
* **Magnetic fields** were considered and not implemented (decision 19).
* **The paper's prose** is a draft produced with AI assistance. Read it against the code and rewrite it in your own words
  before submitting anywhere, and check each venue's rules on AI assistance.

## What I would do next

1. Repeat the experiment in an ephemeris model (Sun + eccentricity) to test whether the k-scaling survives.
2. Add one correction burn and map injection accuracy against burn timing: closer to how missions are flown.
3. Optimise over arrival routes, since that is where the leverage is.
4. Use a realistic injection covariance (burn-direction errors dominate) instead of isotropic Gaussians.
5. Extend to L1 halos and to NRHOs, which need a different capture definition because k approaches 1.
6. Get Richardson (1980) and Howell (1984) through a library and check the coefficients and method against the originals.
