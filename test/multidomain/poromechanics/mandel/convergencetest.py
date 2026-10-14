#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Grid convergence of Mandel's problem.

Runs the executable on successively refined grids and checks the convergence rate of the
discrete L2 errors of pressure and horizontal displacement at the final time.

Usage:
  python3 convergencetest.py <executable> [additional arguments passed to the executable]
"""

import subprocess
import sys
from math import log

LEVELS = [0, 1, 2]
MIN_RATE = 1.9


def collect_errors(executable, args):
    """Run all refinement levels and return the L2 errors of p and u_x per level."""
    errors = {"p": [], "ux": []}
    for level in LEVELS:
        output = subprocess.run(
            [executable, "params.input", *args, "-Grid.Refinement", str(level)],
            check=True, capture_output=True, text=True,
        ).stdout
        line = next(l for l in output.splitlines() if l.startswith("** L2 errors"))
        values = dict(
            entry.split(" = ") for entry in line.split(": ", 1)[1].split(", ")
        )
        errors["p"].append(float(values["L2(p)"]))
        errors["ux"].append(float(values["L2(ux)"]))
    return errors


def rates(errors):
    return [log(errors[i]/errors[i+1])/log(2.0) for i in range(len(errors) - 1)]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    errors = collect_errors(executable, sys.argv[2:])

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
