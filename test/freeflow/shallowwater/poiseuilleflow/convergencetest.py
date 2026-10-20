#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence of the steady Poiseuille flow towards the parabolic velocity profile.

Refines the grid along the channel with a fine grid across the channel, and across the channel
with a fine grid along the channel, and checks the convergence rates of the relative L2 error of
the velocity. The simulation runs until the flow is steady.

With the option --exact-along-channel, the flow is expected to be discretized exactly along the channel
(second-order flux): refined along the channel, the error then stays at the level set by the resolution
across the channel instead of converging.

Usage:
  python3 convergencetest.py <executable> [--exact-along-channel] [additional arguments passed to the executable]
"""

import math
import re
import subprocess
import sys

# (cells along the channel, cells across the channel) and the minimum mean rate of each series
SERIES = {
    "along the channel": ([(10, 160), (20, 160), (40, 160)], 0.85),
    "across the channel": ([(160, 20), (160, 40), (160, 80)], 1.6),
}

STEADY_STATE = ["-TimeLoop.TEnd", "1e6", "-TimeLoop.MaxTimeStepSize", "1e5", "-TimeLoop.PrintoutTimes", "1e6"]

ERRORS = re.compile(r"L1: (\S+) (\S+) (\S+)\n\s*L2: (\S+) (\S+) (\S+)\nNorms of the exact solution \(h, u, v\):\n"
                    r"\s*L1: (\S+) (\S+) (\S+)\n\s*L2: (\S+) (\S+) (\S+)")


def run(executable, cells, args=()):
    """Run the steady channel flow on a grid with the given numbers of cells and return the relative errors"""
    output = subprocess.run(
        [executable, "params.input", *STEADY_STATE, *args,
         "-Problem.Name", f"poiseuilleflow_{cells[0]}x{cells[1]}", "-Grid.Cells", f"{cells[0]} {cells[1]}"],
        check=True, capture_output=True, text=True,
    ).stdout
    values = [float(v) for v in ERRORS.findall(output)[-1]]
    if not all(math.isfinite(v) for v in values):
        raise RuntimeError(f"Non-finite errors for {cells} cells: {values}")
    l2, exactL2 = values[3:6], values[9:12]
    return {"u": l2[1]/exactL2[1]}


def rates(values):
    return [math.log(values[i]/values[i+1])/math.log(2.0) for i in range(len(values) - 1)]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    args = sys.argv[2:]
    exactAlongChannel = "--exact-along-channel" in args
    args = [a for a in args if a != "--exact-along-channel"]

    failed = False
    for name, (grids, minRate) in SERIES.items():
        errors = [run(executable, cells, args)["u"] for cells in grids]
        r = rates(errors)
        print(f"\nrelative L2 error of u, refined {name}:")
        for i, (cells, error) in enumerate(zip(grids, errors)):
            rate = f"{r[i - 1]:.2f}" if i > 0 else "   -"
            print(f"  {cells[0]:>4} x {cells[1]:<4} {error:.4e}  {rate}")
        meanRate = sum(r)/len(r)
        print(f"  mean rate: {meanRate:.2f}")
        if exactAlongChannel and name == "along the channel":
            spread = max(errors)/min(errors) - 1.0
            print(f"  relative spread: {spread:.3f}")
            if spread > 0.1:
                sys.stderr.write(f"Error refined {name} changes by {spread:.3f}, expected a constant error.\n")
                failed = True
        elif meanRate < minRate:
            sys.stderr.write(f"Convergence rate refined {name} ({meanRate:.2f}) is below {minRate}.\n")
            failed = True

    sys.exit(1 if failed else 0)
