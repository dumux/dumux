#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Spatial-convergence test for the coupled Mandel problem.

The executable (test_md_mandel) appends one line per run to two log files
- <flow problem name>.log : L2(p)
- <mech problem name>.log : L2(ux), L2(uy)
This driver runs the executable at several refinement levels, parses the
log files, computes the convergence rates and prints them.
"""

from math import log, isnan, isinf
import os
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write("Please provide the executable name as argument, "
                     "followed by any additional CMD_ARGS (e.g. params.input).\n")
    sys.exit(1)

executableName = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]

# log file names follow the sub-problem names set in params.input
# [OneP.Problem] Name = flow ; [PoroElastic.Problem] Name = mech
flowLog = "flow.log"
mechLog = "mech.log"

# refinement levels passed via -Grid.Refinement <i>; halves h per increment
refinements = [0, 1, 2]

# remove old log files
for f in (flowLog, mechLog):
    subprocess.call(["rm", "-f", f])
    print("Removed old log file ({})!".format(f))

# run the executable at each refinement level
for i in refinements:
    rc = subprocess.call(["./" + executableName] + testargs + ["-Grid.Refinement", str(i)])
    if rc != 0:
        sys.stderr.write("Executable failed at refinement {} (exit code {}).\n".format(i, rc))
        sys.exit(1)

# verify the executable produced the expected log files
for f in (flowLog, mechLog):
    if not os.path.exists(f):
        sys.stderr.write(
            "Expected log file '{}' not found. "
            "Rebuild '{}' so the [ConvergenceTest] output is emitted "
            "(see L2-error block in test_mandel.cc).\n".format(f, executableName))
        sys.exit(1)


def _readErrors(logFile, numErrorColumns):
    """Read a list of error tuples from a [ConvergenceTest] log file."""
    errors = [[] for _ in range(numErrorColumns)]
    with open(logFile, "r") as f:
        for line in f:
            line = line.strip("\n").replace("[ConvergenceTest]", "").split()
            # entries are: name = value name = value ...
            # so the numeric tokens are at positions 2, 5, 8, ...
            for k in range(numErrorColumns):
                errors[k].append(float(line[2 + 3 * k]))
    return errors


def _rate(errPrev, errCurr):
    return (log(errPrev) - log(errCurr)) / log(2.0)


def checkRatesFlow():
    (errorP,) = _readErrors(flowLog, 1)

    rates = []
    with open(flowLog, "w") as f:
        f.write("n\terrorP\t\trateP\n")
        f.write("-" * 50 + "\n")
        for i in range(len(errorP) - 1):
            if isnan(errorP[i]) or isinf(errorP[i]):
                continue
            if errorP[i] < 1e-12 or errorP[i + 1] < 1e-12:
                f.write("error: exact solution!?\n")
                continue
            rateP = _rate(errorP[i], errorP[i + 1])
            f.write("{}\t{:0.4e}\t{:0.4e}\n".format(i, errorP[i], rateP))
            rates.append(rateP)
        i = len(errorP) - 1
        f.write("{}\t{:0.4e}\n".format(i, errorP[i]))

    print("\nComputed convergence rates for pressure:\n")
    subprocess.call(["cat", flowLog])
    return {"p": rates}


def checkRatesMech():
    errorUx, errorUy = _readErrors(mechLog, 2)

    ratesUx, ratesUy = [], []
    with open(mechLog, "w") as f:
        f.write("n\terrorUx\t\trateUx\t\terrorUy\t\trateUy\n")
        f.write("-" * 50 + "\n")
        for i in range(len(errorUx) - 1):
            if isnan(errorUx[i]) or isinf(errorUx[i]):
                continue
            tooSmall = (errorUx[i] < 1e-12 or errorUx[i + 1] < 1e-12) \
                   and (errorUy[i] < 1e-12 or errorUy[i + 1] < 1e-12)
            if tooSmall:
                f.write("error: exact solution!?\n")
                continue
            rUx = _rate(errorUx[i], errorUx[i + 1])
            rUy = _rate(errorUy[i], errorUy[i + 1])
            f.write("{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(
                i, errorUx[i], rUx, errorUy[i], rUy))
            ratesUx.append(rUx)
            ratesUy.append(rUy)
        i = len(errorUx) - 1
        f.write("{}\t{:0.4e}\t\t\t{:0.4e}\n".format(i, errorUx[i], errorUy[i]))

    print("\nComputed convergence rates for displacement:\n")
    subprocess.call(["cat", mechLog])
    return {"ux": ratesUx, "uy": ratesUy}


def reportAll():
    flowRates = checkRatesFlow()
    mechRates = checkRatesMech()

    def mean(xs):
        return float(sum(xs)) / len(xs) if xs else float("nan")

    print()
    for name, rates in (("p",  flowRates["p"]),
                        ("ux", mechRates["ux"]),
                        ("uy", mechRates["uy"])):
        print("Mean rate for {}: {:.3f}".format(name, mean(rates)))


reportAll()
