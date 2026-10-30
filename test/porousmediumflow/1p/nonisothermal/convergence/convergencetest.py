#!/usr/bin/env python3

# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Runs the non-isothermal single-phase convergence test on four refinements and checks
the convergence rate of the temperature in the L2 norm, and that the affine pressure
is represented exactly.
"""

import argparse
import json
import math
import os
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument("executable")
parser.add_argument("name")
parser.add_argument("--expected-rate", type=float, required=True)
args, runtimeArgs = parser.parse_known_args()

if os.path.exists(args.name + ".json"):
    os.remove(args.name + ".json")

for refinement in range(4):
    subprocess.run(
        ["./" + args.executable] + (runtimeArgs or ["params.input"])
        + ["-Problem.Name", args.name, "-Grid.Refinement", str(refinement)],
        check=True,
    )

with open(args.name + ".json") as jsonFile:
    data = json.load(jsonFile)

errors = data["L2errorsTemperature"]
numElements = data["numElements"]
dim = 2
rates = [
    (math.log(errors[i]) - math.log(errors[i + 1]))
    / math.log((numElements[i + 1] / numElements[i]) ** (1.0 / dim))
    for i in range(len(errors) - 1)
]

print("elements  L2 error T   rate")
for i, error in enumerate(errors):
    rate = "{:.3f}".format(rates[i - 1]) if i > 0 else ""
    print("{:8d}  {:.4e}  {}".format(numElements[i], error, rate))
print("L2 error p: " + ", ".join("{:.2e}".format(e) for e in data["L2errorsPressure"]))

failed = False
if min(rates[-2:]) < args.expected_rate - 0.1:
    sys.stderr.write("Temperature converges with rates {} instead of {}\n".format(rates, args.expected_rate))
    failed = True
if max(data["L2errorsPressure"]) > 1e-6:
    sys.stderr.write("The affine pressure is not represented exactly\n")
    failed = True
sys.exit(1 if failed else 0)
