# The paper in plain words

This is the same study as `paper.pdf`, with the technical vocabulary taken out. Use it to understand
the work; use the technical paper when you submit, because examiners expect the proper terms.
Every term in *italics with brackets* is the word the technical paper uses; `docs/GLOSSARY.md` explains each.

## The question

A spacecraft can reach a parking loop near the far side of the Moon by coasting along a special path,
with no braking burn. If the spacecraft starts on that path with a small mistake in position or speed,
does it still arrive? How big a mistake is too big? Does the answer depend on the size of the loop?

## The simplified world

Earth and Moon move on perfect circles and the spacecraft is too light to affect them
(*circular restricted three-body problem*). We watch from a viewpoint that turns with the Moon, so both
bodies stand still. The parking loops are the *L2 halo orbits*; the special coasting paths are the loop's
*stable manifold*.

What this leaves out, and why it matters:

- **The Sun.** Its pull nudges these loops noticeably over a month. Real numbers would differ.
- **The Moon's orbit is slightly oval**, not a circle.
- **The Moon's gravity is lumpy**, which matters most for paths that pass close to it.
- **No corrections.** Real spacecraft fix their course every few days. Here they never do, on purpose,
  so the result is a worst case: how good must you be *before the first correction*?
- **Mistakes are equally likely in every direction.** Real ones are not.

So this is a study of how sensitivity *scales*. It does not design a mission.

## What was done

1. **Built the loops.** A rough formula from 1980 gives a first guess; the computer then adjusts it until
   the path closes exactly, and steps from each loop to the next size. That produced 283 loops.
   Three of them were compared with values published by NASA's Jet Propulsion Laboratory and matched to
   better than 11 digits.
2. **Found the free-ride paths** for each loop, and where they cross the line through the Moon. Those
   crossings form a ring: the target.
3. **Ran the experiment.** For 10 loop sizes, 24 starting points on the ring, 22 sizes of mistake and
   150 random mistakes each: 3.17 million coasting flights. Each was labelled *arrived*, *drifted away* or
   *hit the Moon*.
4. **Used two rules for "arrived"**, because the rule changes the answer:
   - **Rule A:** stays within 5,000 km of the loop for one whole lap.
   - **Rule B:** stays within 124,000 km of the Moon for the same length of time.
5. **Summarised each case by one number:** the mistake size at which half the spacecraft arrive
   (*the 50% tolerance*).
6. **Repeated everything with Einstein's correction to gravity switched on**, using the same random mistakes.

## What was found

| Loop height | Instability number k | 50% tolerance, position (rule A) | 50% tolerance, speed (rule A) |
|---|---|---|---|
| 7,800 km | 588 | about 100 m | 4 mm/s |
| 23,400 km | 462 | about 200 m | 5.5 mm/s |
| 39,000 km | 283 | about 400 m | 9 mm/s |
| 54,600 km | 127 | about 660 m | 16 mm/s |
| 70,100 km | 32 | about 5 km | 69 mm/s |
| 76,000 km | 11 | about 18 km | 185 mm/s |

- **Bigger loops forgive bigger mistakes**, because they are less unstable. Halving the instability number
  raises the position tolerance about 2.3 times.
- **Rule B is 60 to 140 times more forgiving**, and for the two largest loops no 50% point could be found
  at all within the mistakes tested.
- **The starting point on the ring matters more than the loop.** For one loop, the best starting point
  tolerates up to 1,000 times more error than the worst. The fragile ones pass close to the Moon.
- **A simple rule predicts the experiment.** Calculate how many times a mistake grows along the path
  (call it g). Tolerance ≈ 5,000 km ÷ (0.674 × g). This held across four powers of ten.
- **Asking for more laps costs a lot.** Each extra half lap divides the tolerance by about the square
  root of 2k (about 20 for a mid-sized loop).
- **Einstein's correction made no difference.** It changes the force by 4 parts in 100 billion.
  Of 3,168,000 comparisons, 11 came out differently.

## How we know the code is right

Ten independent checks, all passed. Among them: the balance points agree with NASA's values to 13 digits;
a spacecraft placed exactly on L4 does not move for 100 time units; the conserved quantity changes by less
than 1 part in 10 trillion; Mercury's orbit turns by 42.98 arcseconds per century when Einstein's correction
is on and by zero when it is off; GPS clocks come out 38.6 microseconds per day fast. Halving every
accuracy setting changed fewer than 1 in 1,000 individual outcomes and none of the conclusions.

## A mistake that was found and fixed

The first full run failed its own test. Making one internal step ten times smaller changed the answer by
up to four times, which should not happen. The cause: the usual way of picking starting points on the
free-ride tube bunches them unevenly, and the bunching moved when the step changed. Spreading the points
evenly fixed it, and the whole experiment was run again. An odd flat stretch in the first results turned
out to be the same error and disappeared.

## What this does not show

- Real tolerances for a real mission (the Sun and corrections are missing).
- Anything about loops near L1, or the very elongated loops close to the Moon.
- Whether the scaling rule survives in a more realistic model. That is the obvious next experiment.
