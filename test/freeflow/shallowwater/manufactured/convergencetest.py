#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence of the manufactured solution of the shallow water equations.

Refines the grid and checks the convergence rate of the relative L1 errors of water depth and velocity
between the two finest grids.
With --newton-order, the order of convergence of the Newton method is checked as well, from the
maximum relative shifts of a run started from a perturbed exact solution.
With --dt-times-cells C, the time step size is refined with the cell size, dt = C/cells on the unit square.

Usage:
  python3 convergencetest.py <executable> <minimum rate> [--newton-order] [--dt-times-cells C]
                             [additional arguments passed to the executable]
"""

import math
import re
import subprocess
import sys

CELLS = [32, 64, 128]
ERRORS = re.compile(r"L1: (\S+) (\S+) (\S+)\n\s*L2: (\S+) (\S+) (\S+)\nNorms of the exact solution \(h, u, v\):\n"
                    r"\s*L1: (\S+) (\S+) (\S+)\n\s*L2: (\S+) (\S+) (\S+)")
SHIFTS = re.compile(r"Newton iteration\s+\d+ done, maximum relative shift = (\S+)")


def run(executable, cells, args=()):
    output = subprocess.run(
        [executable, "params.input", *args, "-Problem.Name", f"manufactured_{cells}", "-Grid.Cells", f"{cells} {cells}"],
        check=True, capture_output=True, text=True,
    ).stdout
    values = [float(v) for v in ERRORS.findall(output)[-1]]
    if not all(math.isfinite(v) for v in values):
        raise RuntimeError(f"Non-finite errors for {cells} cells: {values}")
    l1, exactL1 = values[0:3], values[6:9]
    return [e/n for e, n in zip(l1, exactL1)], [float(s) for s in SHIFTS.findall(output)]


def rates(values):
    return [math.log(values[i]/values[i+1])/math.log(2.0) for i in range(len(values) - 1)]


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    minRate = float(sys.argv[2])
    args = sys.argv[3:]
    checkNewton = "--newton-order" in args
    args = [a for a in args if a != "--newton-order"]
    dtTimesCells = None
    if "--dt-times-cells" in args:
        i = args.index("--dt-times-cells")
        dtTimesCells = float(args[i + 1])
        del args[i:i + 2]

    def timeStepArgs(cells):
        return ["-TimeLoop.Dt", str(dtTimesCells/cells)] if dtTimesCells is not None else []

    results = [run(executable, cells, [*args, *timeStepArgs(cells)])[0] for cells in CELLS]
    failed = False
    for k, name in enumerate(["relative L1 error of h", "relative L1 error of u", "relative L1 error of v"]):
        errors = [r[k] for r in results]
        r = rates(errors)
        print(f"\n{name}:  cells  error       rate")
        for i, (cells, error) in enumerate(zip(CELLS, errors)):
            rate = f"{r[i - 1]:.2f}" if i > 0 else "   -"
            print(f"       {cells:>5}  {error:.4e}  {rate}")
        if r[-1] < minRate:
            sys.stderr.write(f"Convergence rate of {name} ({r[-1]:.2f}) is below {minRate}.\n")
            failed = True

    if checkNewton:
        _, shifts = run(executable, CELLS[1], [*args, "-Problem.InitialPerturbation", "0.05"])
        print("\nNewton shifts:", " ".join(f"{s:.2e}" for s in shifts))
        orders = [math.log(shifts[i+1])/math.log(shifts[i]) for i in range(len(shifts) - 1)
                  if 1e-10 < shifts[i+1] and shifts[i] < 1e-2]
        print("orders:", " ".join(f"{o:.2f}" for o in orders))
        if not orders or max(orders) < 1.8:
            sys.stderr.write("Newton method does not converge quadratically.\n")
            failed = True

    sys.exit(1 if failed else 0)
