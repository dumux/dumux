#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Trace the wrinkled branch beyond the onset and check that it is a supercritical pitchfork from the
critical stretch found by the eigenvalue analysis: the intercept of the fit of the branch agrees with
the critical stretch, the coefficient c_2 is positive, and the growth stretch grows quadratically
with the amplitude near the onset.

usage: branchtest.py <executable> <parameter file> <amplitudes> <intercept tolerance>
"""

import subprocess
import sys
from math import log


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.stderr.write(__doc__)
        sys.exit(1)

    executable = sys.argv[1] if "/" in sys.argv[1] else "./" + sys.argv[1]
    paramFile, amplitudes, tolerance = sys.argv[2], sys.argv[3], float(sys.argv[4])
    output = subprocess.check_output([executable, paramFile, "-Problem.BranchAmplitudes", amplitudes], text=True)
    print(output)

    values, branch = {}, []
    for line in output.splitlines():
        if line.startswith("Branch amplitude "):
            branch.append((float(line.split()[2].rstrip(":")), float(line.split()[-1])))
        elif ": " in line:
            key, value = line.split(": ", 1)
            values[key] = value

    critical = float(values["Critical growth stretch"])
    deviation = float(values["Branch intercept deviation"])
    c2 = float(values["Branch coefficient c2"])
    (a0, s0), (a1, s1) = branch[0], branch[1]
    exponent = log((s1 - critical)/(s0 - critical))/log(a1/a0)
    print(f"intercept deviation {deviation:.3e}, c2 = {c2:.6f}, exponent {exponent:.4f}")

    if deviation > tolerance or c2 <= 0.0 or abs(exponent - 2.0) > 0.02:
        sys.stderr.write("The branch is not a supercritical pitchfork from the critical stretch. Test failed.\n")
        sys.exit(1)
    sys.exit(0)
