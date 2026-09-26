#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import os
import subprocess
import sys
from math import log


QUANTITIES = {
    "deformation w": "Relative L2-error deformation w: ",
    "in-plane displacement u": "Relative L2-error in-plane displacement u: ",
    "shear gradient potential phi": "Relative L2-error shear gradient potential phi: ",
}


def collect_data(testname, clmax_values, geo="../disk.geo", mesh="disk.msh"):
    """Run the test at each clmax level and return the mesh sizes and the errors it printed."""
    exe = testname if os.path.isabs(testname) else "./" + testname
    hs, errors = [], {name: [] for name in QUANTITIES}
    for level, clmax in enumerate(clmax_values, start=1):
        subprocess.check_call(
            ["gmsh", "-2", "-format", "msh2", "-clmax", str(clmax), geo, "-o", mesh]
        )
        output = subprocess.check_output([exe], text=True)
        for line in output.split("\n"):
            if line.startswith("Max element diameter: "):
                hs.append(float(line.split()[-1]))
            for name, prefix in QUANTITIES.items():
                if line.startswith(prefix):
                    errors[name].append(float(line.split()[-1]))
        active_errors = {name: e for name, e in errors.items() if e}
        if not active_errors or len(hs) != level or any(len(e) != level for e in active_errors.values()):
            raise RuntimeError(f"Missing or inconsistent convergence data from {testname} at clmax={clmax}")
    os.remove(mesh)
    return hs, {name: e for name, e in errors.items() if e}


def compute_rates(hs, errors):
    """Compute convergence rates from actual mesh diameters."""
    if len(hs) < 2 or len(hs) != len(errors):
        raise ValueError("Convergence rates require at least two mesh sizes and one error per mesh")
    return [log(errors[i] / errors[i + 1]) / log(hs[i] / hs[i + 1]) for i in range(len(errors) - 1)]


def print_table(hs, errors, rates):
    """Print a convergence table."""
    print(f"  {'h':>10}  {'L2-error':>12}  {'rate':>6}")
    for i, (h, e) in enumerate(zip(hs, errors)):
        rate_str = f"{rates[i - 1]:.2f}" if i > 0 else "   -"
        print(f"  {h:>10.4f}  {e:>12.4e}  {rate_str:>6}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("testname", help="name of the test executable")
    parser.add_argument("--geo", default="../disk.geo", help="Gmsh geometry file")
    parser.add_argument("--mesh", default="disk.msh", help="mesh file written by Gmsh")
    parser.add_argument("--clmax", type=float, nargs="+", default=[0.25, 0.125, 0.0625])
    args = parser.parse_args()

    testname = args.testname
    hs, errors = collect_data(testname, args.clmax, args.geo, args.mesh)

    failed = False
    for quantity, errs in errors.items():
        rates = compute_rates(hs, errs)
        print(f"\nConvergence rates for {testname} ({quantity}):")
        print_table(hs, errs, rates)

        mean_rate = sum(rates) / len(rates)
        print(f"Mean convergence rate: {mean_rate:.2f}")

        if not (1.8 <= mean_rate <= 2.2):
            sys.stderr.write("*" * 70 + "\n")
            sys.stderr.write(
                f"Convergence rate {mean_rate:.2f} for {quantity} not close enough to 2! Test failed.\n"
            )
            sys.stderr.write("*" * 70 + "\n")
            failed = True

    sys.exit(1 if failed else 0)
