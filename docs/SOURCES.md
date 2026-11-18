# Sources

Rule followed: a source is listed as **opened** only if its page was actually
loaded while this project was built (10 October 2026). Anything quoted without
opening the primary text is marked **not verified at source**, with a note on
how the code is protected against the gap.

## Opened and used

| Source | What was read | Used for |
|---|---|---|
| W. S. Koon, M. W. Lo, J. E. Marsden, S. D. Ross, *Dynamical Systems, the Three-Body Problem and Space Mission Design* (free PDF, v1.3, 2022), https://www.ross.aoe.vt.edu/books/Ross_3BodyProblem_Book_2022.pdf | Chapter 2: Eqs. (2.3.8)-(2.3.10), (2.3.13)-(2.3.14), Table 2.2.1 | Equations of motion, effective potential, energy/Jacobi convention C = -2E |
| NASA/JPL Three-Body Periodic Orbits API v1.0, https://ssd-api.jpl.nasa.gov/periodic_orbits.api and its documentation page | Earth-Moon system block (mu, length and time units, L1, L2) and L2 northern halo rows | mu, units, L-point check, the three published halos of verification check 6 |
| JPL Solar System Dynamics, Astrodynamic Parameters, https://ssd.jpl.nasa.gov/astro_par.html | c, G, au, GM of Sun, Mercury, Earth, Moon, Jupiter system (DE440) | constants.py |
| JPL Solar System Dynamics, Planetary Physical Parameters, https://ssd.jpl.nasa.gov/planets/phys_par.html | radii, rotation periods | constants.py |
| NASA NSSDCA Moon and Mercury fact sheets (D. R. Williams, updated 11 Jan 2024) | Moon radius; Mercury a, e, period | constants.py; Mercury test |
| M. J. Holman et al., "ASSIST: An Ephemeris-Quality Test Particle Integrator", arXiv:2303.16246, Eq. (2) | the full point-mass PPN/EIH acceleration, attributed there to Moyer (2003) | `accel_eih_correction`, `pn_test_particle` |
| N. Ashby, "Relativity in the Global Positioning System", Living Rev. Relativity 6 (2003) 1, section 5, Eqs. (35)-(36) | net fractional offset 4.4647e-10 = 6.9693e-10 - 2.5046e-10 | GPS check |
| Orekit (CS GROUP, Apache-2.0), `org.orekit.orbits.cr3bp.RichardsonExpansion` | all third-order coefficients | cross-check of `richardson.py` (code written independently, compared term by term) |

## Bibliographic record opened, full text not available

| Source | Status |
|---|---|
| D. L. Richardson, "Analytic construction of periodic orbits about the collinear points", Celestial Mechanics 22 (1980) 241-253, doi:10.1007/BF01229511 | Abstract and record opened; article paywalled. **Coefficients not verified at source.** Protected by: Orekit cross-check, and by the fact that the seed is only a starting guess whose corrected result matches JPL. |
| K. C. Howell, "Three-dimensional, periodic, 'halo' orbits", Celestial Mechanics 32 (1984) 53-71, doi:10.1007/BF01358403 | Abstract and record opened. The symmetric-shooting method attributed to it is described from general knowledge. **Not verified at source.** The abstract's statements (families exist for all mass ratios; ranges of stable orbits) are consistent with what this code computes. |
| A. Einstein, L. Infeld, B. Hoffmann, "The Gravitational Equations and the Problem of Motion", Annals of Mathematics 39 (1938) 65-100, doi:10.2307/1968714 | Citation confirmed via the reference list of the Wikipedia article on the EIH equations. Paper itself not opened. The equations used come from Holman et al. (above). |
| V. Szebehely, *Theory of Orbits: The Restricted Problem of Three Bodies*, Academic Press, 1967 | Existence and bibliographic details confirmed via publisher/library listings. **Book not opened.** Cited as the standard reference for the CR3BP and the convention C(L4) = 3; nothing in the code depends on a formula taken only from it. |
| R. S. Park, W. M. Folkner, J. G. Williams, D. H. Boggs, "The JPL Planetary and Lunar Ephemerides DE440 and DE441", AJ 161 (2021) 105 | Citation read from JPL's page; paper not opened. Source of the GM values via JPL's table. |
| E. E. Mamajek et al., IAU 2015 Resolution B3, arXiv:1510.07674 | Abstract opened; it does not state numbers. The nominal solar radius 6.957e8 m is quoted from general knowledge. **Not verified at source.** It enters only the solar-limb deflection check (1.75 arcsec; a 0.1% change in radius changes the result by 0.1%). |

## Not opened at all (general knowledge), with how each is protected

| Item | Protection |
|---|---|
| 1PN conserved energy of the EIH system (Landau & Lifshitz, *Classical Theory of Fields*, section 106) | `test_1pn_energy_is_conserved_better_than_newtonian_energy` |
| Yoshida's 4th-order symplectic coefficients (H. Yoshida, Phys. Lett. A 150 (1990) 262) | `test_yoshida_is_fourth_order` measures the order |
| Schwarzschild geodesic equations, ISCO, photon sphere | checks 8 and `test_schwarzschild_orbits_and_clocks` (weak-field limit recovered) |
| Roche limit coefficients 2^(1/3) (rigid) and 2.44 (fluid); Routh's critical mass ratio | Routh value re-derived in code as (1 - sqrt(23/27))/2 and tested against the eigenvalue computation; Roche coefficients are quoted |
| WGS-84 Earth constants in the GPS check (a1, J2, rotation rate) | result matches Ashby's printed 4.4647e-10 to 5 figures |
| Approximate mass ratios for Sun-Jupiter and Sun-Earth in the interface's system menu | labelled "approximate" in the interface |
| White dwarf, neutron star, black hole presets | labelled "representative" in code and interface |
