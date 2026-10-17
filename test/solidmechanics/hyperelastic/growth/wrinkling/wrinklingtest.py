#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Run the wrinkling onset test on three nested meshes and compare the critical growth stretch and
the number of half wavelengths of the wrinkle with the reference values. The cells of the
coarsest mesh are multiplied by 2 and 4, the grading factors raised to the power 1/2 and 1/4.

usage: wrinklingtest.py <executable> <parameter file> <reference stretch> <reference half wavelengths> <minimum rate>
"""

import subprocess
import sys
from math import copysign, log


def gridParameters(paramFile):
    values, group = {}, ""
    with open(paramFile) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if line.startswith("["):
                group = line.strip("[]")
            elif "=" in line and group == "Grid":
                key, value = (s.strip() for s in line.split("=", 1))
                values[key] = [float(v) for v in value.split()]
    return values


def run(executable, paramFile, grid, factor):
    args = [executable, paramFile]
    for key, value in grid.items():
        if key.startswith("Cells"):
            args += [f"-Grid.{key}", " ".join(str(int(v) * factor) for v in value)]
        elif key.startswith("Grading"):
            args += [f"-Grid.{key}", " ".join(str(copysign(abs(v) ** (1.0 / factor), v)) for v in value)]
    output = subprocess.check_output(args, text=True)
    result = {}
    for line in output.splitlines():
        for key in ["Max element diameter", "Critical growth stretch", "Wrinkle half-wavelengths"]:
            if line.startswith(f"{key}: "):
                result[key] = float(line.split()[-1])
    return result


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    paramFile = sys.argv[2]
    referenceStretch, referenceMode, minRate = float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
    grid = gridParameters(paramFile)
    results = [run(executable, paramFile, grid, factor) for factor in [1, 2, 4]]

    failed = False
    print(f"{'h':>10}{'critical stretch':>20}{'error':>14}{'rate':>8}{'half waves':>12}")
    for i, r in enumerate(results):
        error = abs(r["Critical growth stretch"] - referenceStretch)
        rate = ""
        if i > 0:
            previous = abs(results[i - 1]["Critical growth stretch"] - referenceStretch)
            rate = log(previous / error) / log(results[i - 1]["Max element diameter"] / r["Max element diameter"])
            failed = failed or (i == len(results) - 1 and rate < minRate)
            rate = f"{rate:.2f}"
        mode = int(r["Wrinkle half-wavelengths"])
        failed = failed or mode != referenceMode
        print(f"{r['Max element diameter']:>10.4f}{r['Critical growth stretch']:>20.10f}{error:>14.3e}{rate:>8}{mode:>12}")

    if failed:
        sys.stderr.write(f"Wrong wrinkle mode or convergence rate below {minRate}. Test failed.\n")
        sys.exit(1)
    sys.exit(0)
