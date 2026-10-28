#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Check the clamped square plate against its Legendre-Galerkin reference on refined meshes."""
import os
import re
import subprocess
import sys
from math import log

CELLS = [8, 16, 32, 64, 128]
MIN_ORDER = {"w": 1.8, "M centre": 1.8, "M edge": 0.9, "M edge reaction": 1.8, "phi": 1.8, "psi edge": 0.9}
# the error of psi at the interior point changes sign under refinement, so its value is gated
MAX_ERROR = {"psi interior": 1e-3}


def run(exe, cells):
    out = subprocess.check_output([exe, "params.input", "-Problem.Mapping", "clamped",
                                   "-Grid.Cells", f"{cells} {cells}"], text=True)
    value = lambda key: tuple(map(float, re.search(rf"^{key}: (\S+) reference (\S+)$", out, re.MULTILINE).groups()))
    error = lambda pair: abs(pair[0]/pair[1] - 1.0)
    psi = float(re.search(r"^Largest \|psi\|/\(F L\^2\): (\S+)$", out, re.MULTILINE).group(1))
    return {"w": error(value("Centre deflection")), "M centre": error(value("Centre moment M11")),
            "M edge": error(value("Edge moment M11")),
            "M edge reaction": error(value("Edge moment M11 from the reaction")),
            "phi": error(value("Centre potential phi")), "psi edge": error(value("Edge quarter-point psi")),
            "psi interior": error(value("Interior psi")), "psi": psi}


if __name__ == "__main__":
    exe = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "./test_kirchhoff_love_navier_square")
    results = [run(exe, cells) for cells in CELLS]
    keys = ["w", "M centre", "M edge", "M edge reaction", "phi", "psi edge", "psi interior"]
    print(f"{'cells':>6} " + " ".join(f"{k:>15}" for k in keys) + f" {'max|psi|/(F L^2)':>17}")
    for cells, r in zip(CELLS, results):
        print(f"{cells:>6} " + " ".join(f"{r[k]:15.3e}" for k in keys) + f" {r['psi']:17.4e}")
    failed = False
    for key, bound in MAX_ERROR.items():
        if results[-1][key] > bound:
            sys.stderr.write(f"{key} deviates by {results[-1][key]:.3e} on the finest mesh, more than {bound}\n")
            failed = True
    for key in MIN_ORDER:
        errors = [r[key] for r in results]
        rates = [log(errors[i]/errors[i+1])/log(2.0) for i in range(len(errors) - 1)]
        print(f"  orders in {key}: " + ", ".join(f"{p:.2f}" for p in rates))
        if rates[-1] < MIN_ORDER[key]:
            sys.stderr.write(f"order {rates[-1]:.2f} in {key} below {MIN_ORDER[key]}\n")
            failed = True
    sys.exit(1 if failed else 0)
