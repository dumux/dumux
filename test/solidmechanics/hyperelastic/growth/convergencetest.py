#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Run a growth benchmark on three nested meshes and check the convergence rates of the L2
and H1 errors of the displacement against its exact solution. The cells of the coarsest mesh are
multiplied by 2 and 4 for the finer meshes.

usage: convergencetest.py <executable> <parameter file> <cells of the coarsest mesh> <minimum L2 rate> <minimum H1 rate>
"""

import subprocess
import sys
from math import log

FIELDS = ["L2 error displacement", "H1 error displacement"]


def run(executable, paramFile, cells):
    output = subprocess.check_output(
        [executable, paramFile, "-Grid.Cells", " ".join(str(c) for c in cells)], text=True
    )
    result = {}
    for line in output.splitlines():
        if line.startswith("Max element diameter: "):
            result["h"] = float(line.split()[-1])
        for field in FIELDS:
            if line.startswith(f"{field}: "):
                result[field] = float(line.split()[-1])
    return result


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    paramFile = sys.argv[2]
    coarseCells = [int(c) for c in sys.argv[3].split()]
    minRates = dict(zip(FIELDS, map(float, sys.argv[4:6])))
    results = [run(executable, paramFile, [c * factor for c in coarseCells]) for factor in [1, 2, 4]]

    failed = False
    print(f"{'h':>10}" + "".join(f"{field:>32}" for field in FIELDS))
    for i, r in enumerate(results):
        row = f"{r['h']:>10.4f}"
        for field in FIELDS:
            if i == 0:
                row += f"{r[field]:>24.4e} {'':>7}"
            else:
                rate = log(results[i - 1][field] / r[field]) / log(results[i - 1]["h"] / r["h"])
                failed = failed or rate < minRates[field]
                row += f"{r[field]:>24.4e} ({rate:4.2f})"
        print(row)

    if failed:
        sys.stderr.write(f"A convergence rate is below its minimum {minRates}. Test failed.\n")
        sys.exit(1)
    sys.exit(0)
