#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence for the unsteady 2D Navier-Stokes test of the FVCA8 benchmark.

Runs the executable with central weighting of the advective momentum flux for the viscosity 0.1 on
the coarsest grids of the benchmark and checks the mean convergence orders of the relative
space-time L2 errors of pressure and velocity, as defined by the benchmark, and that the discrete
velocity divergence vanishes.

Usage:
  python3 convergencetest.py <executable>
"""

import csv
import os
import re
import subprocess
import sys
from collections import namedtuple
from math import log, sqrt

CELLS = [4, 8, 16]
VISCOSITY = 0.1
CENTRAL = 0.5
MIN_RATE = 1.8
MAX_DIVERGENCE = 1e-10

Result = namedtuple(
    "Result",
    [
        "name",
        "pressureError",
        "velocityError",
        "newtonIterations",
        "gradientError",
        "divergenceError",
        "nonZeros",
    ],
)


def readCsv(fileName):
    """Read a csv file written by the executable into a list of rows with named columns."""
    with open(fileName, encoding="utf-8") as csvFile:
        rows = list(csv.reader(csvFile, skipinitialspace=True))
    return [dict(zip(rows[0], map(float, row))) for row in rows[1:]]


def timeIntegral(steps, integrand):
    """Integrate a quantity over time with the time step sizes of the given steps."""
    integral, previousTime = 0.0, 0.0
    for step in steps:
        integral += (step["time"] - previousTime) * integrand(step)
        previousTime = step["time"]
    return integral


def spaceTimeError(steps, components):
    """The relative space-time L2 error of the given components."""
    error = timeIntegral(steps, lambda s: sum(s[f"L2Abs({c})"] ** 2 for c in components))
    reference = timeIntegral(
        steps, lambda s: sum((s[f"L2Abs({c})"] / s[f"L2Rel({c})"]) ** 2 for c in components)
    )
    return sqrt(error / reference)


def caseName(cells, viscosity, upwindWeight, averagedVelocity):
    """The name of the output files of one case."""
    return f"angeli_n{cells}_nu{viscosity:g}_w{upwindWeight:g}_a{int(averagedVelocity)}"


def command(executable, cells, viscosity, upwindWeight, averagedVelocity):
    """The command line of one case with the benchmark setup."""
    # the time step size decreases with the square of the cell size, like the expected spatial error
    dt = 0.2 / (viscosity * cells**2)
    parameters = {
        "Problem.Name": caseName(cells, viscosity, upwindWeight, averagedVelocity),
        "Problem.PrintErrors": "true",
        "Problem.FixPressureInOneCell": "true",
        "Problem.ZeroMeanPressureErrors": "true",
        "Problem.InterpolateExactVelocity": str(averagedVelocity).lower(),
        "Vtk.WriteOnlyFinalStep": "true",
        "Component.LiquidKinematicViscosity": str(viscosity),
        "Grid.Cells": f"{cells} {cells}",
        "TimeLoop.DtInitial": str(dt),
        "TimeLoop.MaxTimeStepSize": str(dt),
        "TimeLoop.TEnd": str(1.0 / (10.0 * viscosity)),
        "Flux.UpwindWeight": str(upwindWeight),
    }
    arguments = [executable, "params.input"]
    for key, value in parameters.items():
        arguments += [f"-{key}", value]
    return arguments


def run(executable, cells, viscosity, upwindWeight=CENTRAL, averagedVelocity=True):
    """Run one case and return its errors, Newton iterations and numbers of non-zero entries."""
    name = caseName(cells, viscosity, upwindWeight, averagedVelocity)
    output = subprocess.run(
        command(executable, cells, viscosity, upwindWeight, averagedVelocity),
        check=True,
        stdout=subprocess.PIPE,
        text=True,
        # the direct solver gains little from more threads on these small systems
        env={**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "2")},
    ).stdout

    steps = readCsv(f"{name}_error.csv")
    velocitySteps = readCsv(f"{name}_velocity_error.csv")
    gradientError = timeIntegral(velocitySteps, lambda s: s["gradientError2"])
    gradientReference = timeIntegral(velocitySteps, lambda s: s["gradientReference2"])
    nonZeros = {key: int(value) for key, value in readCsv(f"{name}_matrix.csv")[0].items()}
    newton = len(re.findall(r"Newton iteration\s+\d+ done", output))

    return Result(
        name=name,
        pressureError=spaceTimeError(steps, ["p"]),
        velocityError=spaceTimeError(steps, ["vx", "vy"]),
        newtonIterations=newton / len(steps),
        gradientError=sqrt(gradientError / gradientReference),
        divergenceError=sqrt(timeIntegral(velocitySteps, lambda s: s["divergence2"])),
        nonZeros=nonZeros,
    )


def rates(errors, unknowns):
    """The orders as defined by the benchmark, with respect to the numbers of unknowns."""
    return [
        -2.0 * log(errors[i + 1] / errors[i]) / log(unknowns[i + 1] / unknowns[i])
        for i in range(len(errors) - 1)
    ]


def velocityUnknowns(cells):
    """The numbers of velocity unknowns for the given numbers of cells per direction."""
    return [2 * n * (n + 1) for n in cells]


def pressureUnknowns(cells):
    """The numbers of pressure unknowns for the given numbers of cells per direction."""
    return [n * n for n in cells]


def checkRate(name, errors, unknowns):
    """Print the errors and orders and return whether the mean order is high enough."""
    orders = rates(errors, unknowns)
    print(f"\n{name}:  cells  L2 error    rate")
    for i, (cells, error) in enumerate(zip(CELLS, errors)):
        rate = f"{orders[i - 1]:.2f}" if i > 0 else "   -"
        print(f"       {cells:>5}  {error:.4e}  {rate}")
    meanRate = sum(orders) / len(orders)
    print(f"  mean rate: {meanRate:.2f}")
    if meanRate < MIN_RATE:
        sys.stderr.write(f"Convergence rate of {name} ({meanRate:.2f}) is below {MIN_RATE}.\n")
        return False
    return True


def main():
    """Run the grids, check the orders and the divergence and exit with the result."""
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    results = [run(executable, cells, VISCOSITY) for cells in CELLS]

    passed = checkRate("p", [r.pressureError for r in results], pressureUnknowns(CELLS))
    passed &= checkRate("u", [r.velocityError for r in results], velocityUnknowns(CELLS))

    # the scheme satisfies the discrete mass balance up to round-off
    maxDivergence = max(r.divergenceError for r in results)
    print(f"\nmaximum errdivu: {maxDivergence:.2e}")
    if maxDivergence > MAX_DIVERGENCE:
        sys.stderr.write(f"Velocity divergence ({maxDivergence:.2e}) is above {MAX_DIVERGENCE}.\n")
        passed = False

    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
