#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence of the uniform film flow with viscous bottom friction towards the exact mean velocity.

Refines the grid along the channel and checks the convergence rate of the L2 error of the velocity.

Usage:
  python3 convergencetest.py <executable> [additional arguments passed to the executable]
"""

import math
import re
import subprocess
import sys

CELLS = [20, 40, 80]
MIN_RATE = 0.9

ERRORS = re.compile(r"L1: (\S+) (\S+) (\S+)\n\s*L2: (\S+) (\S+) (\S+)")


def run(executable, cells, args=()):
    """Run the film flow on a grid with the given number of cells along the channel and return the errors"""
    output = subprocess.run(
        [executable, "params.input", *args, "-Problem.Name", f"filmflow_{cells}", "-Grid.Cells", f"{cells} 2"],
        check=True, capture_output=True, text=True,
    ).stdout
    errors = [float(e) for e in ERRORS.findall(output)[-1]]
    if not all(math.isfinite(e) for e in errors):
        raise RuntimeError(f"Non-finite errors for {cells} cells: {errors}")
    return {"L2(h)": errors[3], "L2(u)": errors[4]}


def rates(values):
    return [math.log(values[i]/values[i+1])/math.log(2.0) for i in range(len(values) - 1)]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    errors = [run(executable, cells, sys.argv[2:])["L2(u)"] for cells in CELLS]
    r = rates(errors)

    print("\nL2(u):  cells  error       rate")
    for i, (cells, error) in enumerate(zip(CELLS, errors)):
        rate = f"{r[i - 1]:.2f}" if i > 0 else "   -"
        print(f"        {cells:>5}  {error:.4e}  {rate}")
    meanRate = sum(r)/len(r)
    print(f"  mean rate: {meanRate:.2f}")

    if meanRate < MIN_RATE:
        sys.stderr.write(f"Convergence rate ({meanRate:.2f}) is below {MIN_RATE}.\n")
        sys.exit(1)
    sys.exit(0)
