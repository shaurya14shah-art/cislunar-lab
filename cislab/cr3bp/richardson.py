"""Richardson's third-order analytic approximation of halo orbits.

Source: D. L. Richardson, "Analytic construction of periodic orbits about the
collinear points", Celestial Mechanics 22 (1980) 241-253, doi:10.1007/BF01229511.

Honesty note.  The journal article is paywalled; its bibliographic record and
abstract were opened, but the coefficient formulas below could not be read from
the article itself.  They were therefore cross-checked, term by term, against
an independent open-source implementation (Orekit, org.orekit.orbits.cr3bp.
RichardsonExpansion, Apache-2.0), and - more importantly - they are *tested*:
the seed they produce is fed to the differential corrector, which converges to
orbits that match NASA/JPL's published catalogue (see tests/test_halo.py).

Idea of the method (Lindstedt-Poincare)
---------------------------------------
Move the origin to the Lagrange point and measure lengths in units of gamma
(its distance from the nearer primary).  Expand the gravitational potential in
Legendre polynomials; the coefficients are

    c_n = (1/gamma^3) [ (+-1)^n mu + (-1)^n (1-mu) gamma^(n+1) / (1 -+ gamma)^(n+1) ]

(upper signs L1, lower signs L2).  The linear solution has *different* in-plane
(lambda) and out-of-plane (sqrt(c2)) frequencies, so it does not close.  A halo
orbit appears when the amplitudes are large enough for the nonlinear terms to
pull the two frequencies together; that requirement is the amplitude constraint

    l1 Ax^2 + l2 Az^2 + Delta = 0,     Delta = lambda^2 - c2.

The solution is then built as a Fourier series in tau1 = lambda*nu*t up to third
order in the amplitudes, with the frequency correction nu = 1 + s1 Ax^2 + s2 Az^2.
"""
from __future__ import annotations
import numpy as np
from .lagrange import collinear_gamma


class Richardson:
    def __init__(self, mu, point="L2"):
        if point not in ("L1", "L2"):
            raise ValueError("Richardson expansion is implemented for L1 and L2")
        self.mu, self.point = mu, point
        g = self.gamma = collinear_gamma(mu, point)

        def c(n):
            if point == "L1":
                return (mu + (-1) ** n * (1 - mu) * g ** (n + 1) / (1 - g) ** (n + 1)) / g**3
            return ((-1) ** n * mu + (-1) ** n * (1 - mu) * g ** (n + 1) / (1 + g) ** (n + 1)) / g**3

        c2, c3, c4 = c(2), c(3), c(4)
        self.c2, self.c3, self.c4 = c2, c3, c4
        # lam: in-plane linear frequency (root of lam^4 + (c2-2) lam^2 - (c2-1)(1+2c2) = 0)
        lam = self.lam = np.sqrt(0.5 * (2 - c2 + np.sqrt(9 * c2 * c2 - 8 * c2)))
        # k: ratio of y to x amplitude of the linear in-plane ellipse
        k = self.k = (lam * lam + 1 + 2 * c2) / (2 * lam)
        d1 = 3 * lam**2 / k * (k * (6 * lam**2 - 1) - 2 * lam)
        d2 = 8 * lam**2 / k * (k * (11 * lam**2 - 1) - 2 * lam)
        a21 = 3 * c3 * (k * k - 2) / (4 * (1 + 2 * c2))
        a22 = 3 * c3 / (4 * (1 + 2 * c2))
        a23 = -3 * c3 * lam / (4 * k * d1) * (3 * k**3 * lam - 6 * k * (k - lam) + 4)
        a24 = -3 * c3 * lam / (4 * k * d1) * (2 + 3 * k * lam)
        b21 = -3 * c3 * lam / (2 * d1) * (3 * k * lam - 4)
        b22 = 3 * c3 * lam / d1
        d21 = -c3 / (2 * lam**2)
        a31 = (-9 * lam / (4 * d2) * (4 * c3 * (k * a23 - b21) + k * c4 * (4 + k * k))
               + (9 * lam**2 + 1 - c2) / (2 * d2) * (2 * c3 * (2 * a23 - k * b21) + c4 * (2 + 3 * k * k)))
        a32 = (-9 * lam / (4 * d2) * (4 * c3 * (k * a24 - b22) + k * c4)
               - 3 * (9 * lam**2 + 1 - c2) / (2 * d2) * (c3 * (k * b22 + d21 - 2 * a24) - c4))
        b31 = (3 / (8 * d2) * 8 * lam * (3 * c3 * (k * b21 - 2 * a23) - c4 * (2 + 3 * k * k))
               + 3 / (8 * d2) * (9 * lam**2 + 1 + 2 * c2) * (4 * c3 * (k * a23 - b21) + k * c4 * (4 + k * k)))
        b32 = (9 * lam / d2 * (c3 * (k * b22 + d21 - 2 * a24) - c4)
               + 3 / (8 * d2) * (9 * lam**2 + 1 + 2 * c2) * (4 * c3 * (k * a24 - b22) + k * c4))
        d31 = 3 / (64 * lam**2) * (4 * c3 * a24 + c4)
        d32 = 3 / (64 * lam**2) * (4 * c3 * (a23 - d21) + c4 * (4 + k * k))
        den = 2 * lam * (lam * (1 + k * k) - 2 * k)
        s1 = (1.5 * c3 * (2 * a21 * (k * k - 2) - a23 * (k * k + 2) - 2 * k * b21)
              - 0.375 * c4 * (3 * k**4 - 8 * k * k + 8)) / den
        s2 = (1.5 * c3 * (2 * a22 * (k * k - 2) + a24 * (k * k + 2) + 2 * k * b22 + 5 * d21)
              + 0.375 * c4 * (12 - k * k)) / den
        l1 = -1.5 * c3 * (2 * a21 + a23 + 5 * d21) - 0.375 * c4 * (12 - k * k) + 2 * lam**2 * s1
        l2 = 1.5 * c3 * (a24 - 2 * a22) + 1.125 * c4 + 2 * lam**2 * s2
        self.coef = dict(a21=a21, a22=a22, a23=a23, a24=a24, a31=a31, a32=a32, b21=b21, b22=b22,
                         b31=b31, b32=b32, d21=d21, d31=d31, d32=d32, s1=s1, s2=s2, l1=l1, l2=l2)
        self.delta = lam * lam - c2

    # -- amplitudes are in units of gamma inside this class ------------------
    def ax_from_az(self, az):
        """In-plane amplitude forced by the amplitude constraint."""
        q = (-self.delta - self.coef["l2"] * az * az) / self.coef["l1"]
        if q < 0:
            raise ValueError("no halo orbit of this Az exists in the third-order theory")
        return float(np.sqrt(q))

    @property
    def min_ax(self):
        """Smallest in-plane amplitude at which halos bifurcate from the planar family."""
        return float(np.sqrt(-self.delta / self.coef["l1"]))

    def period(self, az):
        ax = self.ax_from_az(az)
        nu = 1 + self.coef["s1"] * ax * ax + self.coef["s2"] * az * az
        return float(2 * np.pi / (self.lam * nu))

    def state(self, az_nd, t=0.0, north=True):
        """Approximate barycentric rotating-frame state on the halo.

        ``az_nd`` is the out-of-plane amplitude in CR3BP length units (NOT gamma
        units); multiply by the Earth-Moon distance to get kilometres.
        """
        q = self.coef
        g, lam, k = self.gamma, self.lam, self.k
        az = az_nd / g
        ax = self.ax_from_az(az)
        nu = 1 + q["s1"] * ax * ax + q["s2"] * az * az
        w = lam * nu
        t1 = w * t
        dm = 1.0 if north else -1.0          # Richardson's delta_m = 2 - m, m = 1 or 3
        c1, s1_, c2_, s2_, c3_, s3_ = (np.cos(t1), np.sin(t1), np.cos(2 * t1), np.sin(2 * t1),
                                       np.cos(3 * t1), np.sin(3 * t1))
        x = (q["a21"] * ax * ax + q["a22"] * az * az - ax * c1
             + (q["a23"] * ax * ax - q["a24"] * az * az) * c2_
             + (q["a31"] * ax**3 - q["a32"] * ax * az * az) * c3_)
        y = (k * ax * s1_ + (q["b21"] * ax * ax - q["b22"] * az * az) * s2_
             + (q["b31"] * ax**3 - q["b32"] * ax * az * az) * s3_)
        z = dm * (az * c1 + q["d21"] * ax * az * (c2_ - 3) + (q["d32"] * az * ax * ax - q["d31"] * az**3) * c3_)
        vx = w * (ax * s1_ - 2 * (q["a23"] * ax * ax - q["a24"] * az * az) * s2_
                  - 3 * (q["a31"] * ax**3 - q["a32"] * ax * az * az) * s3_)
        vy = w * (k * ax * c1 + 2 * (q["b21"] * ax * ax - q["b22"] * az * az) * c2_
                  + 3 * (q["b31"] * ax**3 - q["b32"] * ax * az * az) * c3_)
        vz = dm * w * (-az * s1_ - 2 * q["d21"] * ax * az * s2_
                       - 3 * (q["d32"] * az * ax * ax - q["d31"] * az**3) * s3_)
        x_l = 1 - self.mu + (g if self.point == "L2" else -g)
        return np.array([x * g + x_l, y * g, z * g, vx * g, vy * g, vz * g])
