#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Grid convergence of axisymmetric Stokes flow with a manufactured solution (nonzero radial
velocity), for the symmetrized and the unsymmetrized viscous stress and for a domain on the
axis (0 < r < 1) and an annulus (0.5 < r < 1.5). Usage:
convergencetest.py <executable> <minimum L2 rate> [further arguments for the executable].
The L2 error of the velocity must decrease at least with the given rate between the two
finest grids.
"""

import math
import os
import subprocess
import sys

executable, minimumRate = sys.argv[1], float(sys.argv[2])
extraArgs = sys.argv[3:]
inputFile = os.path.join(os.path.dirname(os.path.abspath(executable)), "params.input")
cells = [8, 16, 32, 64]


def run(numCells, args):
    output = subprocess.run(
        [executable, inputFile, "-Grid.Cells", f"{numCells} {numCells}", *args, *extraArgs],
        check=True, capture_output=True, text=True,
    ).stdout
    tokens = next(l for l in output.splitlines() if l.startswith("[ConvergenceTest]")).split()[1:]
    return {tokens[i]: float(tokens[i + 1]) for i in range(0, len(tokens), 2)}


cases = {
    "axis, symmetrized": ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "false"],
    "axis, unsymmetrized": ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"],
    "annulus, symmetrized": ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "false",
                             "-Grid.LowerLeft", "0 0.5", "-Grid.UpperRight", "1 1.5"],
    "annulus, unsymmetrized": ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true",
                               "-Grid.LowerLeft", "0 0.5", "-Grid.UpperRight", "1 1.5"],
}

failed = False
for name, args in cases.items():
    results = [run(n, args) for n in cells]
    print(f"{name}:")
    print(f"{'cells':>8}{'L2 error':>16}{'rate':>8}{'H1 error':>16}{'rate':>8}")
    for i, r in enumerate(results):
        rates = [math.log(results[i - 1][k]/r[k])/math.log(2.0) if i > 0 else float("nan") for k in ["errorL2", "errorH1"]]
        print(f"{int(r['cells']):8d}{r['errorL2']:16.6e}{rates[0]:8.3f}{r['errorH1']:16.6e}{rates[1]:8.3f}")
    if rates[0] < minimumRate:
        failed = True

sys.exit(1 if failed else 0)
