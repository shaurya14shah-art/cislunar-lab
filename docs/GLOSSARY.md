# Glossary

Every technical term used in the project, in one line each.

| Term | Meaning |
|---|---|
| **50% tolerance, σ50** | The mistake size at which half the spacecraft still arrive. |
| **Amplitude, Az** | The size of a halo: its greatest height above or below the Moon's orbital plane. |
| **Ballistic capture** | Arriving in an orbit by coasting, with no braking burn. |
| **Bootstrap range** | An estimate of how much a result might change if the experiment were repeated, made by reshuffling the trials already run. |
| **Capture curve** | This project's name for the ring where the free-ride tube crosses the line through the Moon. |
| **Centrifugal force** | The outward push you feel in a turning viewpoint. |
| **Circular restricted three-body problem (CR3BP)** | The simplified model used throughout: Earth and Moon move on perfect circles, and the spacecraft is too light to disturb them. |
| **Common random numbers** | Using the same random draws in two experiments so that any difference comes from the thing being changed, not from luck. |
| **Continuation** | Stepping from one loop to the next size, using each as the starting guess for its neighbour. |
| **Coriolis force** | The sideways push on anything moving in a turning viewpoint. |
| **Differential correction** | Repeatedly adjusting a first guess until the path closes exactly. Newton's method applied to an orbit. |
| **Eigenvalue** | A number that says how fast a small disturbance grows or how fast it oscillates. A positive real one means runaway growth. |
| **Einstein–Infeld–Hoffmann equations** | The 1PN equations of motion for several bodies. |
| **Error growth factor, g** | How many times larger a starting mistake has become by the end of the test. Tolerance ≈ 5,000 km ÷ (0.674 g). |
| **Exponential growth** | Growing by the same factor in every equal stretch of time. |
| **Floquet multiplier** | The exact factor by which a mistake grows in one lap. The largest one is about 2k. |
| **Halo orbit** | A tilted closed loop near L1 or L2. From Earth it appears to circle the Moon like a halo. |
| **Hill radius** | The distance from the Moon within which the Moon's gravity dominates the Earth's: about 62,000 km. |
| **Hill region** | The set of places a spacecraft with a given C is able to reach. Water on the chart. |
| **Injection** | The moment the spacecraft is placed on its coasting path. Injection error is the mistake made then. |
| **Integrator** | The procedure that advances a path forward in small time steps. DOP853 is the high-accuracy one used here. |
| **ISCO** | Innermost stable circular orbit: the closest a particle can circle a black hole without falling in. Three Schwarzschild radii. |
| **Jacobi constant, C** | The one quantity that stays fixed along any path in this model. It acts as energy turned upside down: larger C, less energy. |
| **Lagrange point (L1–L5)** | A balance point where a spacecraft at rest stays at rest in the turning view. |
| **Logarithmic scale** | An axis where each step multiplies instead of adds: 1, 10, 100, 1000. |
| **Lyapunov orbit** | The flat cousin of a halo orbit, lying in the Moon's orbital plane. |
| **Mass ratio, μ** | The Moon's share of the total Earth-plus-Moon mass: 0.01215, about 1.2%. The only number the model needs. |
| **Monodromy matrix** | The state transition matrix for exactly one lap of a loop. |
| **Monte Carlo** | Answering a question by running many trials with random inputs and counting outcomes. |
| **Near-rectilinear halo orbit (NRHO)** | A very elongated, almost stable halo that passes close to the Moon. Planned home of the lunar Gateway station. |
| **Non-dimensional units** | Measuring length in Earth–Moon distances and time so that one month is 2π. It makes the equations free of awkward constants. |
| **Orbit family** | The whole continuous range of loops of different sizes. |
| **Periodic orbit** | A path that closes on itself and repeats. |
| **Poincaré section** | Marking only the points where paths cross a chosen surface, to turn a tangle of curves into a readable picture. |
| **Post-Newtonian (1PN)** | Newton's gravity plus the first and largest correction from Einstein's general relativity. |
| **Power law** | A relationship of the form y = constant × x to some power. It appears as a straight line when both axes are logarithmic. |
| **Richardson approximation** | A formula from 1980 that gives a rough first guess for a halo orbit. |
| **Roche limit** | The distance inside which a large body's tides would pull a small one apart. |
| **Rotating frame** | A viewpoint that turns with the Moon, so Earth and Moon appear fixed. |
| **Schwarzschild radius** | The size to which a mass must be squeezed to become a black hole. For Earth, 9 mm. |
| **Stability index, k** | The instability number. A mistake grows about 2k times per lap. k = 1 means stable. |
| **Stable manifold** | The tube of free-ride paths that glide onto an unstable loop. |
| **Standard deviation, σ (sigma)** | The typical size of a random mistake. About two-thirds of mistakes are smaller than one σ. |
| **State** | Position and velocity together: six numbers that fix what the spacecraft does next. |
| **State transition matrix (STM)** | A table that answers: if I change the start slightly, how does the end change? It is how error growth is calculated without random trials. |
| **Symplectic integrator** | A stepping method that never lets energy drift away over long runs. |
| **Tidal force** | The stretching that comes from gravity being stronger on the near side of an object than the far side. |
| **Time dilation** | Clocks tick slower when moving fast or sitting deep in a gravity well. |
| **Unstable / stable** | Unstable: a small miss grows (hilltop). Stable: a small miss stays small (bowl). |
| **Unstable manifold** | The tube of paths that leave the loop. |
| **Zero-velocity curve / surface** | The coastline: the edge of the region a spacecraft can reach. It would have to be at rest to touch it. |
