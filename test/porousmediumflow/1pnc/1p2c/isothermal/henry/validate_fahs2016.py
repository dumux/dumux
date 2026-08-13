#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Runs the Henry-problem simulation given by --command, then validates its output
against the digitized semianalytical isochlor positions of Fahs et al. (2016, WRR,
doi:10.1002/2016WR019288), Appendix D, Table D1 ("Test Case 1": the classic, purely
diffusive Henry [1964] problem).

This is a physics validation against an independent, externally published reference
solution (a different numerical method entirely: Fourier-Galerkin semianalytical, not
just a re-run of our own code) -- not a byte-for-byte regression test. The default
--tolerance was empirically tightened (this repository's convention, see e.g. the
lockexchange test) after observing max relative errors of 0.0186 (Test Case 1) and
0.0223 (Test Case 2) at the current 240x80 grid/1 d TEnd (see params.input,
params_case2.input); the default leaves roughly 2x headroom above the worse of the two.
"""

import argparse
import shlex
import subprocess
import sys
import csv

import numpy as np
import meshio


def readReferenceTable(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(row for row in f if not row.startswith("#")):
            rows.append(
                {
                    "Z": float(row["Z"]),
                    10: float(row["X10"]),
                    50: float(row["X50"]),
                    90: float(row["X90"]),
                }
            )
    return rows


def extractIsochlorX(points, values, z, levels, atol=1e-6):
    """Given the scalar field `values` at `points` (an unstructured but grid-aligned
    point cloud), extract the x-position of each dimensionless-concentration level in
    `levels` along the horizontal line at the given z, via linear interpolation."""
    mask = np.isclose(points[:, 1], z, atol=atol)
    if not np.any(mask):
        return {level: None for level in levels}

    x = points[mask, 0]
    c = np.ravel(values[mask])
    order = np.argsort(x)
    x, c = x[order], c[order]

    result = {}
    for level in levels:
        target = level / 100.0
        if target < c.min() or target > c.max():
            result[level] = None
        else:
            # np.interp expects the independent variable (here c) increasing;
            # the concentration profile along a horizontal line is expected to be
            # monotonically increasing from the freshwater to the seawater side
            result[level] = float(np.interp(target, c, x))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--command", required=True, help="The simulation executable and arguments, as a single string")
    parser.add_argument("--vtu", required=True, help="Path to the resulting VTU file to validate")
    parser.add_argument("--reference", required=True, help="Path to the digitized reference CSV (Z,X10,X50,X90)")
    parser.add_argument("--field", default="X^solute_liq", help="VTU point data field name for the salt mass fraction")
    parser.add_argument("--seawater-fraction", type=float, default=0.035, help="Mass fraction at c=1 (seawater) in this simulation's FluidSystem")
    parser.add_argument("--tolerance", type=float, default=0.05, help="Maximum allowed relative error in isochlor x-position (empirically set -- see module docstring)")
    args = parser.parse_args()

    try:
        result = subprocess.call(shlex.split(args.command))
    except OSError:
        print(f"OSError: could not run command: {args.command}")
        sys.exit(1)
    if result:
        sys.exit(result)

    reference = readReferenceTable(args.reference)

    mesh = meshio.read(args.vtu)
    points = mesh.points
    if args.field not in mesh.point_data:
        print(f"Field '{args.field}' not found in {args.vtu}.")
        print(f"Available point data fields: {list(mesh.point_data)}")
        sys.exit(1)
    dimensionlessC = mesh.point_data[args.field] / args.seawater_fraction

    levels = [10, 50, 90]
    header = f"{'Z':>6} {'level':>6} {'X_sim':>10} {'X_ref':>10} {'abs err':>10} {'rel err':>10}"
    print(header)
    print("-" * len(header))

    maxRelError = 0.0
    numMissing = 0
    for row in reference:
        z = row["Z"]
        simX = extractIsochlorX(points, dimensionlessC, z, levels)
        for level in levels:
            xRef = row[level]
            xSim = simX[level]
            if xSim is None:
                print(f"{z:6.3f} {level:5d}% {'--':>10} {xRef:10.3f} {'--':>10} {'--':>10}")
                numMissing += 1
                continue
            absErr = abs(xSim - xRef)
            relErr = absErr / abs(xRef)
            maxRelError = max(maxRelError, relErr)
            print(f"{z:6.3f} {level:5d}% {xSim:10.3f} {xRef:10.3f} {absErr:10.4f} {relErr:10.4f}")

    print("-" * len(header))
    print(f"Max relative error: {maxRelError:.4f} (tolerance: {args.tolerance:.4f})")
    if numMissing:
        print(f"WARNING: {numMissing} isochlor position(s) could not be extracted (level not spanned at that Z).")

    if numMissing or maxRelError > args.tolerance:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
