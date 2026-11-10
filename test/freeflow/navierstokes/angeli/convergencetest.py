#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence for the unsteady 2D Navier-Stokes test of the FVCA8 benchmark.

Runs the executable with central weighting of the advective momentum flux for the viscosity 0.1 on
the coarsest grids of the benchmark and checks the mean convergence orders of the relative
space-time L2 errors of pressure and velocity, as defined by the benchmark.

Usage:
  python3 convergencetest.py <executable> [additional arguments passed to the executable]
"""

import csv
import os
import re
import subprocess
import sys
from math import log, sqrt

CELLS = [4, 8, 16]
VISCOSITY = 0.1
CENTRAL = 0.5
MIN_RATE = 1.8
MAX_DIVERGENCE = 1e-10


def run(executable, cells, viscosity, upwindWeight=CENTRAL, averagedVelocity=True, args=()):
    """Run one case and return its name, the relative space-time L2 errors of pressure and velocity and the Newton iterations per time step."""
    # the time step size decreases with the square of the cell size, like the expected spatial error
    dt = 0.2/(viscosity*cells**2)
    name = f"angeli_n{cells}_nu{viscosity:g}_w{upwindWeight:g}_a{int(averagedVelocity)}"
    result = subprocess.run(
        [executable, "params.input", *args,
         "-Problem.Name", name,
         "-Problem.PrintErrors", "true",
         "-Problem.FixPressureInOneCell", "true",
         "-Problem.ZeroMeanPressureErrors", "true",
         "-Problem.InterpolateExactVelocity", str(averagedVelocity).lower(),
         "-Vtk.WriteOnlyFinalStep", "true",
         "-Component.LiquidKinematicViscosity", str(viscosity),
         "-Grid.Cells", f"{cells} {cells}",
         "-TimeLoop.DtInitial", str(dt),
         "-TimeLoop.MaxTimeStepSize", str(dt),
         "-TimeLoop.TEnd", str(1.0/(10.0*viscosity)),
         "-Flux.UpwindWeight", str(upwindWeight)],
        check=True, stdout=subprocess.PIPE, text=True,
        # the direct solver gains little from more threads on these small systems
        env={**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "2")},
    )
    with open(f"{name}_error.csv") as errorFile:
        rows = list(csv.reader(errorFile, skipinitialspace=True))
    steps = [dict(zip(rows[0], map(float, row))) for row in rows[1:]]

    def spaceTimeError(components):
        error, reference, previousTime = 0.0, 0.0, 0.0
        for step in steps:
            dtStep = step["time"] - previousTime
            previousTime = step["time"]
            for c in components:
                error += dtStep*step[f"L2Abs({c})"]**2
                reference += dtStep*(step[f"L2Abs({c})"]/step[f"L2Rel({c})"])**2
        return sqrt(error/reference)

    with open(f"{name}_velocity_error.csv") as errorFile:
        rows = list(csv.reader(errorFile, skipinitialspace=True))
    velocitySteps = [dict(zip(rows[0], map(float, row))) for row in rows[1:]]
    gradientError, gradientReference, divergence, previousTime = 0.0, 0.0, 0.0, 0.0
    for step in velocitySteps:
        dtStep = step["time"] - previousTime
        previousTime = step["time"]
        gradientError += dtStep*step["gradientError2"]
        gradientReference += dtStep*step["gradientReference2"]
        divergence += dtStep*step["divergence2"]

    with open(f"{name}_matrix.csv") as matrixFile:
        rows = list(csv.reader(matrixFile, skipinitialspace=True))
    nonZeros = dict(zip(rows[0], map(int, rows[1])))

    newton = len(re.findall(r"Newton iteration\s+\d+ done", result.stdout))
    return (name, spaceTimeError(["p"]), spaceTimeError(["vx", "vy"]), newton/len(steps),
            sqrt(gradientError/gradientReference), sqrt(divergence), nonZeros)


def rates(errors, unknowns):
    """The orders as defined by the benchmark, with respect to the numbers of unknowns."""
    return [-2.0*log(errors[i+1]/errors[i])/log(unknowns[i+1]/unknowns[i]) for i in range(len(errors) - 1)]


def velocityUnknowns(cells):
    return [2*n*(n + 1) for n in cells]


def pressureUnknowns(cells):
    return [n*n for n in cells]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    results = [run(executable, cells, VISCOSITY, args=sys.argv[2:]) for cells in CELLS]
    errors = {"p": ([r[1] for r in results], pressureUnknowns(CELLS)),
              "u": ([r[2] for r in results], velocityUnknowns(CELLS))}

    failed = False
    for name, (values, unknowns) in errors.items():
        r = rates(values, unknowns)
        print(f"\n{name}:  cells  L2 error    rate")
        for i, (cells, error) in enumerate(zip(CELLS, values)):
            rate = f"{r[i - 1]:.2f}" if i > 0 else "   -"
            print(f"       {cells:>5}  {error:.4e}  {rate}")
        meanRate = sum(r)/len(r)
        print(f"  mean rate: {meanRate:.2f}")
        if meanRate < MIN_RATE:
            sys.stderr.write(f"Convergence rate of {name} ({meanRate:.2f}) is below {MIN_RATE}.\n")
            failed = True

    # the scheme satisfies the discrete mass balance up to round-off
    maxDivergence = max(r[5] for r in results)
    print(f"\nmaximum errdivu: {maxDivergence:.2e}")
    if maxDivergence > MAX_DIVERGENCE:
        sys.stderr.write(f"Velocity divergence ({maxDivergence:.2e}) is above {MAX_DIVERGENCE}.\n")
        failed = True

    sys.exit(1 if failed else 0)
