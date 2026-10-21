#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence of the wet dam break towards the exact solution of the Riemann problem.

Runs the executable on successively refined grids with a time step size proportional to the
cell size and checks the convergence rate of the relative L1 errors of water depth and velocity
at the final time.

Usage:
  python3 convergencetest.py <executable> [additional arguments passed to the executable]
"""

import math
import re
import subprocess
import sys

CELLS = [100, 200, 400]
MIN_RATE = 0.7

# time step size on the coarsest grid of params.input (100 cells)
TIME_STEP_SIZE_100 = 0.005

ERRORS = re.compile(r"L1: (\S+) (\S+) (\S+)\n\s*L2: (\S+) (\S+) (\S+)\nNorms of the exact solution \(h, u, v\):\n"
                    r"\s*L1: (\S+) (\S+) (\S+)\n\s*L2: (\S+) (\S+) (\S+)")


def run(executable, cells, args=(), name="dambreak"):
    """Run the dam break on a grid with the given number of cells along the channel and return the relative errors"""
    dt = TIME_STEP_SIZE_100*100/cells
    output = subprocess.run(
        [executable, "params.input", *args, "-Problem.Name", f"{name}_{cells}",
         "-Grid.Cells0", str(cells), "-Grid.Cells1", "1",
         "-TimeLoop.MaxTimeStepSize", str(dt), "-TimeLoop.DtInitial", str(dt)],
        check=True, capture_output=True, text=True,
    ).stdout
    values = [float(v) for v in ERRORS.findall(output)[-1]]
    if not all(math.isfinite(v) for v in values):
        raise RuntimeError(f"Non-finite errors for {cells} cells: {values}")
    l1, exactL1 = values[0:3], values[6:9]
    return {"h": l1[0]/exactL1[0], "u": l1[1]/exactL1[1]}


def rates(values):
    return [math.log(values[i]/values[i+1])/math.log(2.0) for i in range(len(values) - 1)]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    results = [run(executable, cells, sys.argv[2:]) for cells in CELLS]

    failed = False
    for name in ["h", "u"]:
        errors = [r[name] for r in results]
        r = rates(errors)
        print(f"\nrelative L1 error of {name}:  cells  error       rate")
        for i, (cells, error) in enumerate(zip(CELLS, errors)):
            rate = f"{r[i - 1]:.2f}" if i > 0 else "   -"
            print(f"                            {cells:>5}  {error:.4e}  {rate}")
        meanRate = sum(r)/len(r)
        print(f"  mean rate: {meanRate:.2f}")
        if meanRate < MIN_RATE:
            sys.stderr.write(f"Convergence rate of {name} ({meanRate:.2f}) is below {MIN_RATE}.\n")
            failed = True

    sys.exit(1 if failed else 0)
