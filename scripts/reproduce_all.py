"""ONE command that regenerates everything the paper reports.

    python scripts/reproduce_all.py            # reuses data/mc if it exists
    python scripts/reproduce_all.py --fresh    # deletes data/mc and reruns the Monte Carlo (~2-3 h on 2 cores)

Steps: Monte Carlo -> analysis -> figures -> numbers for LaTeX -> verification report -> paper PDF.
"""
import os, shutil, subprocess, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
py = sys.executable


def run(*cmd, cwd=ROOT, check=True):
    print("\n>>>", " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=cwd, check=check)


if "--fresh" in sys.argv:
    shutil.rmtree(os.path.join(ROOT, "data", "mc"), ignore_errors=True)
run(py, "scripts/run_experiment.py")
run(py, "scripts/analyze_experiment.py")
if "--fresh" in sys.argv or not os.path.exists(os.path.join(ROOT, "data", "route_study.json")):
    run(py, "scripts/route_study.py")
run(py, "paper/make_figures.py")
run(py, "paper/make_numbers.py")
run(py, "scripts/verification_report.py")
if shutil.which("pdflatex"):
    for _ in range(2):
        run("pdflatex", "-interaction=nonstopmode", "-halt-on-error", "paper.tex", cwd=os.path.join(ROOT, "paper"))
    print("\nPaper: paper/paper.pdf")
else:
    print("\npdflatex not found: figures and numbers are ready; compile paper/paper.tex with any LaTeX distribution.")
