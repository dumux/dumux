#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence of the Taylor-Couette flow.

Runs the executable on successively refined grids and checks the convergence rate of the
discrete L2 errors of pressure and velocity.

Usage:
  python3 convergencetest.py <executable> [additional arguments passed to the executable]
"""

import csv
import subprocess
import sys
from math import log, sqrt

LEVELS = [0, 1, 2]
MIN_RATE = 1.9


def run(executable, level, args=()):
    """Run one refinement level and return its name and the L2 errors of pressure and velocity."""
    name = f"taylorcouette_level{level}"
    subprocess.run(
        [executable, "params.input", *args, "-Problem.Name", name, "-Grid.Refinement", str(level)],
        check=True, stdout=subprocess.DEVNULL,
    )
    with open(f"{name}_error.csv") as errorFile:
        rows = list(csv.reader(errorFile, skipinitialspace=True))
    errors = dict(zip(rows[0], map(float, rows[-1])))
    return name, errors["L2Abs(p)"], sqrt(errors["L2Abs(vx)"]**2 + errors["L2Abs(vy)"]**2)


def rates(errors):
    return [log(errors[i]/errors[i+1])/log(2.0) for i in range(len(errors) - 1)]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    results = [run(executable, level, sys.argv[2:]) for level in LEVELS]
    errors = {"p": [r[1] for r in results], "u": [r[2] for r in results]}

    failed = False
    for name, values in errors.items():
        r = rates(values)
        print(f"\n{name}:  level  L2 error    rate")
        for level, error in zip(LEVELS, values):
            rate = f"{r[level - 1]:.2f}" if level > 0 else "   -"
            print(f"       {level:>5}  {error:.4e}  {rate}")
        meanRate = sum(r)/len(r)
        print(f"  mean rate: {meanRate:.2f}")
        if meanRate < MIN_RATE:
            sys.stderr.write(f"Convergence rate of {name} ({meanRate:.2f}) is below {MIN_RATE}.\n")
            failed = True

    sys.exit(1 if failed else 0)
