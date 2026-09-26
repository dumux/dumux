#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
"""Check simply supported square plate convergence against Navier's series."""
import os
import re
import subprocess
import sys
from math import isfinite, log

CELLS = [8, 16, 32, 64, 128]
MIN_ORDER = {"w": 1.8, "M12": 1.3}


def parse_results(out):
    def values(pattern, count):
        matches = list(re.finditer(pattern, out, re.MULTILINE))
        if len(matches) != count:
            raise ValueError(f"Expected {count} result(s) matching {pattern!r}, got {len(matches)}")
        result = [tuple(map(float, match.groups())) for match in matches]
        if not all(isfinite(value) for row in result for value in row):
            raise ValueError(f"Nonfinite result matching {pattern!r}")
        return result

    def relative_error(value, reference):
        if reference == 0.0:
            raise ValueError("Expected a nonzero Navier reference value")
        return abs(value/reference - 1.0)

    w = values(r"^Centre deflection: (\S+) reference (\S+)$", 1)[0]
    m = values(r"^Centre moment M11: (\S+) reference (\S+)$", 1)[0]
    corners = values(r"^Corner \(([^,]+),([^)]+)\) M12: (\S+) reference (\S+)$", 4)
    if {(x, y) for x, y, _, _ in corners} != {(0, 0), (1, 0), (1, 1), (0, 1)}:
        raise ValueError("Expected one result at each corner of the unit square")
    maxError = values(r"^Maximum deflection error: (\S+)$", 1)[0][0]
    result = {
        "w": relative_error(*w),
        "wmax": maxError,
        "M11": relative_error(*m),
        "M12": max(relative_error(c, r) for _, _, c, r in corners),
    }
    if not all(isfinite(error) and error >= 0.0 for error in result.values()):
        raise ValueError("Expected finite, nonnegative errors")
    return result


def run(exe, mapping, cells):
    out = subprocess.check_output([exe, "params.input", "-Problem.Mapping", mapping,
                                   "-Grid.Cells", f"{cells} {cells}"], text=True)
    return parse_results(out)


def orders(errors):
    if len(errors) < 2 or not all(isfinite(error) and error > 0.0 for error in errors):
        raise ValueError("Convergence rates require at least two finite, positive errors")
    return [(log(errors[i]) - log(errors[i+1]))/log(2.0) for i in range(len(errors) - 1)]


if __name__ == "__main__":
    exe = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "./test_kirchhoff_love_navier_square")

    failed = False
    for mapping in ["tangential", "traction"]:
        results = [run(exe, mapping, cells) for cells in CELLS]
        print(f"\n{mapping} mapping, relative deviations from Navier's series")
        print(f"{'cells':>6} {'w centre':>10} {'max |w|':>10} {'M11 centre':>11} {'M12 corners':>12}")
        for cells, r in zip(CELLS, results):
            print(f"{cells:>6} {r['w']:10.3e} {r['wmax']:10.3e} {r['M11']:11.3e} {r['M12']:12.3e}")
        for key in ["w", "M12"]:
            rates = orders([r[key] for r in results])
            print(f"  orders in {key}: " + ", ".join(f"{p:.2f}" for p in rates))
            if mapping == "tangential" and rates[-1] < MIN_ORDER[key]:
                sys.stderr.write(f"tangential mapping: order {rates[-1]:.2f} in {key} below {MIN_ORDER[key]}\n")
                failed = True

    sys.exit(1 if failed else 0)
