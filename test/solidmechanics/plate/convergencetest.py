#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import sys
import subprocess
from math import log


def collect_data(testname, clmax_values):
    """Run the test at each clmax level and return (hs, errors)."""
    exe = testname if os.path.isabs(testname) else "./" + testname
    hs, errors = [], []
    for clmax in clmax_values:
        subprocess.call(
            ["gmsh", "-2", "-format", "msh2", "-clmax", str(clmax), "../disk.geo", "-o", "disk.msh"]
        )
        output = subprocess.check_output([exe], text=True)
        for line in output.split("\n"):
            if line.startswith("Max element diameter: "):
                hs.append(float(line.split()[-1]))
            if line.startswith("Relative L2-error deformation w: "):
                errors.append(float(line.split()[-1]))
    os.remove("disk.msh")
    return hs, errors


def compute_rates(hs, errors):
    """Compute convergence rates from actual mesh diameters."""
    return [log(errors[i] / errors[i + 1]) / log(hs[i] / hs[i + 1]) for i in range(len(errors) - 1)]


def print_table(hs, errors, rates):
    """Print a convergence table."""
    print(f"  {'h':>10}  {'L2-error':>12}  {'rate':>6}")
    for i, (h, e) in enumerate(zip(hs, errors)):
        rate_str = f"{rates[i - 1]:.2f}" if i > 0 else "   -"
        print(f"  {h:>10.4f}  {e:>12.4e}  {rate_str:>6}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("Please provide the test executable name as argument\n")
        sys.exit(1)

    testname = str(sys.argv[1])
    hs, errors = collect_data(testname, [0.25, 0.125, 0.0625])
    rates = compute_rates(hs, errors)

    print(f"\nConvergence rates for {testname}:")
    print_table(hs, errors, rates)

    mean_rate = sum(rates) / len(rates)
    print(f"Mean convergence rate: {mean_rate:.2f}")

    if not (1.8 <= mean_rate <= 2.2):
        sys.stderr.write("*" * 70 + "\n")
        sys.stderr.write(f"Convergence rate {mean_rate:.2f} not close enough to 2! Test failed.\n")
        sys.stderr.write("*" * 70 + "\n")
        sys.exit(1)

    sys.exit(0)
