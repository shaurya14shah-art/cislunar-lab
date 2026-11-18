# Cislunar Lab

A verified simulator for the circular restricted three-body problem (CR3BP), a
general N-body sandbox with first-order relativity, and an interactive 3D
"navigation chart" front end. Built to support the research paper in `paper/`:

> For an Earth-Moon L2 halo orbit, how large an injection error can a spacecraft
> absorb before stable-manifold capture fails, and how does that tolerance
> scale with orbit amplitude?

## Quick start

Needs Python 3.12 or 3.13 for the exact pinned versions (built and verified on 3.13).
On Python 3.10 or 3.11 use `requirements-min.txt` instead and run the tests to confirm.

```bash
cd cislunar-lab
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python -m cislab.server          # then open http://127.0.0.1:8000
```

The first start compiles the numerical kernels (about 20 seconds); after that
it is instant. Everything runs on your computer; no internet is needed.

```bash
pytest tests                              # 50 automated tests, about a minute
python scripts/verification_report.py     # the ten mandatory checks -> docs/VERIFICATION.md
python scripts/reproduce_all.py           # ONE command: experiment, analysis, figures, report, paper PDF
```

`reproduce_all.py` reuses the Monte Carlo results shipped in `data/mc/`. Add
`--fresh` to rerun the experiment from nothing (two to three hours on two cores).

## Using the interface

The scene fills the window. There are two sheets, switched in the title block.

**Three-body chart.** Water is where a spacecraft of the chosen energy may go;
land is forbidden. Drag the *Jacobi constant* slider to raise or lower the sea.
Click a Lagrange point for its eigenvalues and a live instability experiment.
Click *Halo orbit* to change its size, grow manifolds, fly it, run a live
injection-error Monte Carlo, compare Newtonian with 1PN, or open a Poincare
section. The chart key (bottom left) toggles layers; each row shows the real
line style.

**N-body sandbox.** Pick a scenario or build your own. Click a body: its
properties open beside it and can be edited; the regime meter shows how
relativistic it is. Magenta overlays are relativistic: clock dials, light rays,
the Newtonian-against-1PN ghost trails. *Close-up* opens exact strong-field
orbits and light rays around a single mass.

The ruler at the bottom plays, reverses, changes speed and scrubs time, and
always shows the conserved-quantity drift. *Save picture* writes a PNG, *Save
data* a CSV. Space bar plays/pauses; Esc closes notes.

## Using the physics core directly

```python
from cislab.constants import MU_EARTH_MOON as mu
from cislab.cr3bp import lagrange_points, linear_frequencies
from cislab.cr3bp.families import build_family, member_with_z
from cislab.cr3bp.manifolds import manifold, branch_towards_secondary

print(lagrange_points(mu)["L2"], linear_frequencies(mu, "L2"))
family = build_family(mu)                        # ~280 L2 halo orbits in a few seconds
halo = member_with_z(family, 0.10)               # far-side z = 0.10 (39,000 km)
print(halo.period, halo.jacobi, halo.stability_index)
tube = manifold(halo, "stable", n_phases=40, eps=1e-6,
                sign=branch_towards_secondary(halo, "stable"), section_x=1 - mu)

from cislab.nbody.system import two_body
s = two_body("Sun", "Mercury", 5.79e10)
hist = s.run(3.0e7, method="dop853", pn=True)    # 1PN (Einstein-Infeld-Hoffmann)
print(abs(hist["energy_drift"]).max())
```

## Layout

```
cislab/                 physics core (no rendering code anywhere in here)
  constants.py          sourced constants and body presets
  cr3bp/                Mode A: equations, Lagrange points, Hill regions, Richardson seed,
                        differential correction, continuation, manifolds, Poincare, Monte Carlo
  nbody/                Mode B: Newtonian + 1PN force kernels, Yoshida and adaptive integrators, frames
  relativity/           regime meter, EIH experiments, Schwarzschild geodesics, clocks, light bending
  fields/               potential, field lines, tides, Hill spheres, Roche limits
  server/               FastAPI + WebSocket layer (contains no physics)
web/                    browser front end (Three.js, vendored; works offline)
tests/                  pytest suite; test_verification.py holds the ten mandatory checks
scripts/                run_experiment, analyze_experiment, verification_report, reproduce_all
paper/                  paper.tex, make_figures.py, make_numbers.py, figures/
data/mc/                raw Monte Carlo outcome tables; data/mc_summary.json the analysed results
docs/                   HOW_IT_WORKS, VERIFICATION, DECISIONS, SOURCES, SUMMARY
```

## Start here if the vocabulary is new

* `docs/PLAIN_GUIDE.md`: the whole project in everyday language, with all eight graphs explained.
* `docs/GLOSSARY.md`: every technical term in one line.
* `paper/PLAIN_LANGUAGE_VERSION.md`: the paper without the jargon.
* `web/playground.html`: a simple flat simulation. Double-click it; it needs no server. (Also at `/playground.html` when the server is running.)

## Read next

* `docs/HOW_IT_WORKS.md`: the physics and the code in plain language.
* `docs/VERIFICATION.md`: each check with its measured value.
* `docs/SOURCES.md`: every reference, and which ones were actually opened.
* `docs/DECISIONS.md`: why things are the way they are.
* `docs/SUMMARY.md`: what is done, what is verified, what is not, what next.

## Honest scope

This is a research and teaching tool for an idealised model. It does not design
real missions: the CR3BP has no Sun, no lunar eccentricity, no lumpy gravity and
no correction burns. The relativity layer is first post-Newtonian for N bodies
and exact only for a single non-rotating mass; the interface says which regime
each tool is valid in.

## Licences

Project code: MIT (see `LICENSE`). Vendored: Three.js (MIT, `web/vendor/THREE-LICENSE`);
Alegreya and Alegreya Sans typefaces (SIL Open Font License, `web/fonts/OFL-LICENSE`).
