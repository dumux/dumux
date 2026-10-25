#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Run a convergence test on successively refined grids and check the rates of the L2 error
of the pressure and of the H1 seminorm error of its gradient against the largest element length.

Usage: convergencetest.py <executable> <name> --l2-rate <r> --h1-rate <r>
                          [--levels <n>] [runtime arguments]
"""

import json
import math
import os
import subprocess
import sys

args = sys.argv[1:]


def popOption(name, default):
    if name in args:
        i = args.index(name)
        value = args[i + 1]
        del args[i : i + 2]
        return value
    return default


expectedL2 = float(popOption("--l2-rate", 2.0))
expectedH1 = float(popOption("--h1-rate", 1.0))
levels = int(popOption("--levels", 4))
executable, name = args[0], args[1]
runtimeArgs = args[2:] if len(args) > 2 else ["params.input"]

if os.path.exists(name + ".json"):
    os.remove(name + ".json")

for level in range(levels):
    subprocess.run(
        ["./" + executable] + runtimeArgs + ["-Problem.Name", name, "-Grid.Refinement", str(level)],
        check=True,
    )

with open(name + ".json") as jsonFile:
    data = json.load(jsonFile)

h = data["h"]
failed = False
for key, expected in (("L2errors", expectedL2), ("H1errors", expectedH1)):
    errors = data[key]
    rates = [math.log(errors[i] / errors[i + 1]) / math.log(h[i] / h[i + 1]) for i in range(len(h) - 1)]
    print(f"{key}: " + " ".join(f"{e:.4e}" for e in errors))
    print(f"{key} rates: " + " ".join(f"{r:.3f}" for r in rates))
    if abs(rates[-1] - expected) > 0.1:
        sys.stderr.write(f"{key}: finest rate {rates[-1]:.3f} differs from the expected {expected}\n")
        failed = True

sys.exit(1 if failed else 0)
