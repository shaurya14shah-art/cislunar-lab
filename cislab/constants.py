"""Physical constants and body presets, each with the source it was read from.

Every value below was copied from a page that was opened while this project was
built (October 2026).  ``SOURCES`` maps a short key to the reference, and each
preset names the key(s) it relies on.  Where a number is *representative* rather
than a measurement of one specific object (the neutron star and the black hole)
that is said explicitly, because pretending otherwise would be dishonest.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math

SOURCES = {
    "JPL-ASTRO": "JPL Solar System Dynamics, 'Astrodynamic Parameters', "
                 "https://ssd.jpl.nasa.gov/astro_par.html (c, G: CODATA 2018; GM values: "
                 "Park et al. 2021, AJ 161:105, DE440).",
    "JPL-PHYS": "JPL Solar System Dynamics, 'Planetary Physical Parameters', "
                "https://ssd.jpl.nasa.gov/planets/phys_par.html (radii: Archinal et al. 2018).",
    "NSSDC-MOON": "NASA NSSDCA Moon Fact Sheet (D. R. Williams, updated 11 Jan 2024), "
                  "https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html",
    "NSSDC-MERCURY": "NASA NSSDCA Mercury Fact Sheet (updated 11 Jan 2024), "
                     "https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html",
    "JPL-3BP": "NASA/JPL Three-Body Periodic Orbits API v1.0, "
               "https://ssd-api.jpl.nasa.gov/periodic_orbits.api (Earth-Moon system block).",
    "IAU-2015-B3": "IAU 2015 Resolution B3 nominal solar radius 6.957e8 m "
                   "(Mamajek et al., arXiv:1510.07674). NOTE: only the abstract was opened; "
                   "the numerical value is quoted from general knowledge and is marked unverified.",
    "ASHBY-2003": "N. Ashby, 'Relativity in the Global Positioning System', "
                  "Living Rev. Relativity 6 (2003) 1, Eqs. (35)-(36).",
    "REPRESENTATIVE": "Representative textbook-scale value, not a measurement of one object.",
}

# --- fundamental constants (JPL-ASTRO) -------------------------------------
C = 299_792_458.0            # speed of light, m/s (exact, SI definition)
G = 6.67430e-11              # Newtonian constant, m^3 kg^-1 s^-2 (CODATA 2018)
AU = 149_597_870_700.0       # astronomical unit, m (IAU 2012 B1, exact)
DAY = 86_400.0               # s
JULIAN_CENTURY = 36_525.0 * DAY
ARCSEC = math.pi / (180.0 * 3600.0)   # radians in one arcsecond

# --- gravitational parameters GM, m^3/s^2 (JPL-ASTRO, DE440) ----------------
GM_SUN = 1.32712440041279419e20
GM_MERCURY = 22031.868551e9
GM_EARTH = 398600.435507e9
GM_MOON = 4902.800118e9
GM_JUPITER_SYSTEM = 126712764.1e9

# --- Earth-Moon CR3BP scaling (JPL-3BP) -------------------------------------
MU_EARTH_MOON = 1.215058560962404e-2      # m_Moon / (m_Earth + m_Moon)
LUNIT_EM_KM = 389703.264829278            # length unit: Earth-Moon distance, km
TUNIT_EM_S = 382981.289129055             # time unit: 1/(mean motion), s
VUNIT_EM_KMS = LUNIT_EM_KM / TUNIT_EM_S   # ~1.0176 km/s
R_MOON_KM = 1737.4                        # NSSDC-MOON volumetric mean radius
R_EARTH_KM = 6371.0084                    # JPL-PHYS mean radius
C_ND_EM = (C / 1000.0) / VUNIT_EM_KMS     # speed of light in Earth-Moon units (~2.95e5)


@dataclass
class Body:
    """One gravitating body for the N-body sandbox (SI units throughout).

    ``spin`` is the angular velocity vector in rad/s.  It is carried for display
    and for the rotational-flattening / frame-dragging *estimates*; it does not
    enter the point-mass equations of motion (and the docs say so).
    """
    name: str
    gm: float                                  # G*M, m^3/s^2 (better known than M itself)
    radius: float                              # m
    pos: tuple = (0.0, 0.0, 0.0)               # m
    vel: tuple = (0.0, 0.0, 0.0)               # m/s
    spin: tuple = (0.0, 0.0, 0.0)              # rad/s
    source: str = "REPRESENTATIVE"
    note: str = ""

    @property
    def mass(self) -> float:
        return self.gm / G

    @property
    def schwarzschild_radius(self) -> float:
        """r_s = 2GM/c^2: the radius inside which nothing escapes."""
        return 2.0 * self.gm / C**2

    def to_dict(self):
        d = asdict(self)
        d["mass"] = self.mass
        d["rs"] = self.schwarzschild_radius
        return d


def _spin(period_s):
    return (0.0, 0.0, 2 * math.pi / period_s)


PRESETS = {
    "Sun": Body("Sun", GM_SUN, 6.957e8, spin=_spin(25.38 * DAY), source="JPL-ASTRO; IAU-2015-B3",
                note="Radius is the IAU nominal value (unverified here); rotation ~25.4 d is "
                     "the usual sidereal equatorial figure and is representative."),
    "Mercury": Body("Mercury", GM_MERCURY, 2439.4e3, spin=_spin(58.6462 * DAY),
                    source="JPL-ASTRO; JPL-PHYS"),
    "Earth": Body("Earth", GM_EARTH, 6371.0084e3, spin=_spin(0.99726968 * DAY),
                  source="JPL-ASTRO; JPL-PHYS"),
    "Moon": Body("Moon", GM_MOON, 1737.4e3, spin=_spin(655.720 * 3600.0),
                 source="JPL-ASTRO; NSSDC-MOON"),
    "Jupiter": Body("Jupiter", GM_JUPITER_SYSTEM, 69911e3, spin=_spin(0.41354 * DAY),
                    source="JPL-ASTRO; JPL-PHYS",
                    note="GM is for the whole Jupiter system (planet + moons), as tabulated by JPL."),
    "White dwarf": Body("White dwarf", 1.0 * GM_SUN, 5.8e6, source="REPRESENTATIVE",
                        note="1 solar mass, ~0.0084 solar radii: typical of a Sirius-B-like "
                             "white dwarf. Representative, not a specific measurement."),
    "Neutron star": Body("Neutron star", 1.4 * GM_SUN, 12e3, spin=_spin(0.03), source="REPRESENTATIVE",
                         note="Canonical 1.4 solar masses, 12 km radius, 30 ms spin. "
                              "1PN is NOT trustworthy near this body."),
    "Black hole (10 Msun)": Body("Black hole (10 Msun)", 10.0 * GM_SUN, 2 * 10.0 * GM_SUN / C**2,
                                 source="REPRESENTATIVE",
                                 note="Stellar-mass, non-rotating. 'Radius' is the Schwarzschild "
                                      "radius 2GM/c^2 (the event horizon)."),
}

# Mercury orbit elements used by the perihelion-precession test (NSSDC-MERCURY)
MERCURY_A = 57.909e9          # semi-major axis, m
MERCURY_E = 0.2056            # eccentricity
MERCURY_PERIOD_D = 87.969     # sidereal period, days

# GPS numbers used by the clock test (ASHBY-2003)
GPS_A = 26_561_750.0          # GPS semi-major axis, m (Ashby 2003)
EARTH_GEOID_POTENTIAL_OVER_C2 = 6.969290134e-10   # L_G, defining constant (IAU 2000 B1.9)
