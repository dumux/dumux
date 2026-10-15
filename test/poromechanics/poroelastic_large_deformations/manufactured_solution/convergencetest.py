#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Run the manufactured solution on three nested meshes and check that the L2 errors of
displacement, solid pressure and fluid pressure decrease at second order.
"""

import sys
import subprocess
from math import log

FIELDS = ["displacement", "solid pressure", "fluid pressure"]
MIN_RATE = 1.8


def run(executable, cells):
    output = subprocess.check_output(
        [executable, "params.input", "-Grid.Cells", f"{cells} {cells}"], text=True
    )
    result = {}
    for line in output.splitlines():
        if line.startswith("Max element diameter: "):
            result["h"] = float(line.split()[-1])
        for field in FIELDS:
            if line.startswith(f"L2 error {field}: "):
                result[field] = float(line.split()[-1])
    return result


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("Please provide the test executable as argument\n")
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    results = [run(executable, cells) for cells in [8, 16, 32]]

    failed = False
    print(f"{'h':>10}" + "".join(f"{field:>26}" for field in FIELDS))
    for i, r in enumerate(results):
        row = f"{r['h']:>10.4f}"
        for field in FIELDS:
            if i == 0:
                row += f"{r[field]:>18.4e} {'':>7}"
            else:
                rate = log(results[i - 1][field] / r[field]) / log(results[i - 1]["h"] / r["h"])
                failed = failed or rate < MIN_RATE
                row += f"{r[field]:>18.4e} ({rate:4.2f})"
        print(row)

    if failed:
        sys.stderr.write(f"A convergence rate is below {MIN_RATE}. Test failed.\n")
        sys.exit(1)
    sys.exit(0)
