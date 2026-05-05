#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import subprocess
import sys
from math import log

if len(sys.argv) < 2:
    sys.stderr.write("Please provide the test executable name as argument\n")
    sys.exit(1)

testname = str(sys.argv[1])

errors = []
hs = [0.25, 0.125, 0.0625]
for h in hs:
    subprocess.call(
        ["gmsh", "-2", "-format", "msh2", "-clmax", str(h), "disk.geo", "-o", "disk.msh"]
    )
    output = subprocess.check_output(["./" + testname], text=True)
    for line in output.split("\n"):
        if line.startswith("Relative L2-error deformation w: "):
            errors.append(float(line.split()[-1]))

# cleanup generated mesh file
os.remove("disk.msh")

rates = [(log(errors[i]) - log(errors[i + 1])) / log(2) for i in range(len(errors) - 1)]

print("\nConvergence rates for {}:".format(testname))
for i, (h, e) in enumerate(zip(hs, errors)):
    rate_str = "{:.2f}".format(rates[i - 1]) if i > 0 else "-"
    print("  h={:.4f}  L2-error={:.4e}  rate={}".format(h, e, rate_str))

mean_rate = sum(rates) / len(rates)
print("Mean convergence rate: {:.2f}".format(mean_rate))

if not (1.8 <= mean_rate <= 2.2):
    sys.stderr.write("*" * 70 + "\n")
    sys.stderr.write(
        "Convergence rate {:.2f} not close enough to 2! Test failed.\n".format(mean_rate)
    )
    sys.stderr.write("*" * 70 + "\n")
    sys.exit(1)

sys.exit(0)
