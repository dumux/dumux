#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Runs the heat pipe test on successively refined grids and checks the convergence
rates of the discrete L2 errors with respect to the solution of Udell (1985).

usage: convergencetest.py <executable> [input file] [runtime arguments]
"""

import os
import sys
import subprocess
from math import log

FIELDS = ["saturation", "gas pressure", "liquid pressure", "temperature"]


def collect_data(testname, refinements, input_file="params.input", extra_args=()):
    """Run the test for each grid refinement and return the mesh sizes and the errors per field."""
    exe = testname if os.path.isabs(testname) else "./" + testname
    name = os.path.basename(testname)
    error_file = name + "-errors.csv"
    if os.path.exists(error_file):
        os.remove(error_file)
    for refinement in refinements:
        subprocess.run(
            [exe, input_file, "-Problem.Name", name, "-Grid.Refinement", str(refinement), *extra_args],
            check=True,
            stdout=subprocess.DEVNULL,
        )
    with open(error_file) as f:
        rows = [[float(v) for v in line.split(",")] for line in f if line.strip()]
    hs = [row[0] for row in rows]
    errors = {field: [row[i + 1] for row in rows] for i, field in enumerate(FIELDS)}
    return hs, errors


def compute_rates(hs, errors):
    """Compute convergence rates from the mesh sizes."""
    return [log(errors[i] / errors[i + 1]) / log(hs[i] / hs[i + 1]) for i in range(len(errors) - 1)]


def print_table(hs, errors, rates):
    """Print a convergence table with one error and rate column per field."""
    print(f"  {'h':>10}" + "".join(f"  {field:>16}  {'rate':>6}" for field in FIELDS))
    for i, h in enumerate(hs):
        line = f"  {h:>10.5f}"
        for field in FIELDS:
            rate = f"{rates[field][i - 1]:.2f}" if i > 0 else "-"
            line += f"  {errors[field][i]:>16.4e}  {rate:>6}"
        print(line)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    testname = sys.argv[1]
    input_file = sys.argv[2] if len(sys.argv) > 2 else "params.input"
    hs, errors = collect_data(testname, refinements=range(4), input_file=input_file, extra_args=sys.argv[3:])
    rates = {field: compute_rates(hs, errors[field]) for field in FIELDS}

    print(f"\nConvergence rates for {testname}:")
    print_table(hs, errors, rates)

    # full upwinding of the mobilities makes the scheme first order
    expected_order = 1
    tol = 0.2
    failed = False
    for field in FIELDS:
        mean_rate = sum(rates[field]) / len(rates[field])
        print(f"Mean convergence rate {field}: {mean_rate:.2f}")
        if not (expected_order - tol <= mean_rate <= expected_order + tol):
            sys.stderr.write(f"Convergence rate {mean_rate:.2f} of the {field} error not close enough to {expected_order}!\n")
            failed = True

    sys.exit(1 if failed else 0)
