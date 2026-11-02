#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

from math import isinf, isnan, log
import os
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write("Please provide a single argument <testname> to the script\n")
    sys.exit(1)

testname = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]


def parseRates(testargs):
    """Remove the -Rates "ErrorType minRate maxRate" entries from the arguments and return them"""
    if "-Rates" not in testargs:
        sys.stderr.write("Please provide the expected rates with -Rates \"ErrorType minRate maxRate\" ...\n")
        sys.exit(1)

    idx = testargs.index("-Rates")
    end = idx + 1
    while end < len(testargs) and not testargs[end].startswith("-"):
        end += 1

    rates = {}
    for spec in testargs[idx + 1:end]:
        parts = spec.split()
        if len(parts) != 3:
            sys.stderr.write("Each -Rates entry must be 'ErrorType minRate maxRate', got: '{}'\n".format(spec))
            sys.exit(1)
        rates[parts[0]] = (float(parts[1]), float(parts[2]))
    return rates, testargs[:idx] + testargs[end:]


rateSpecs, testargs = parseRates(testargs)

refinements = [0, 1, 2]

# the error types that can be checked: the csv file and the column of the error
errorColumns = {
    "L2(v)": (testname + "_errors_velocity.csv", 2),
    "H1(v)": (testname + "_errors_velocity.csv", 3),
    "L2(p)": (testname + "_errors_pressure.csv", 2),
    "H1(p)": (testname + "_errors_pressure.csv", 3),
}
unknownRates = set(rateSpecs) - set(errorColumns)
if unknownRates:
    sys.stderr.write("Unknown error types in -Rates: {}\n".format(", ".join(sorted(unknownRates))))
    sys.exit(1)
quantities = [(label, *errorColumns[label], bounds) for label, bounds in rateSpecs.items()]

for fileName in [testname + ".log", testname + "_errors_velocity.csv", testname + "_errors_pressure.csv"]:
    if os.path.exists(fileName):
        os.remove(fileName)

for refinement in refinements:
    subprocess.run(
        ["./" + testname] + testargs + [
            "-Grid.Refinement", str(refinement),
            "-Problem.Name", testname,
            "-Problem.PrintErrors", "false",
            "-Problem.PrintConvergenceTestFile", "true",
        ],
        check=True,
    )


def readCSV(fileName):
    with open(fileName, "r") as csvFile:
        return [
            [float(value.strip()) for value in line.split(",")]
            for line in csvFile
            if line.strip()
        ]


def convergenceRates(rows, errorColumn):
    rates = []
    for i in range(len(rows) - 1):
        h0, h1 = rows[i][1], rows[i + 1][1]
        e0, e1 = rows[i][errorColumn], rows[i + 1][errorColumn]
        if isnan(e0) or isinf(e0) or e0 < 1e-15 or e1 < 1e-15:
            continue
        rates.append(log(e0/e1)/log(h0/h1))
    return rates


data = {}
for _, fileName, _, _ in quantities:
    if fileName not in data:
        data[fileName] = readCSV(fileName)

failed = False
with open(testname + ".log", "w") as logFile:
    logFile.write("quantity\trates\texpected\n")
    logFile.write("-" * 50 + "\n")
    for label, fileName, errorColumn, (minRate, maxRate) in quantities:
        rates = convergenceRates(data[fileName], errorColumn)
        ratesString = ", ".join("{:.4e}".format(rate) for rate in rates)
        logFile.write("{}\t{}\t[{}, {}]\n".format(label, ratesString, minRate, maxRate))
        if len(rates) == 0 or not minRate <= rates[-1] <= maxRate:
            failed = True

print("\nComputed the following convergence rates for {}:\n".format(testname))
subprocess.run(["cat", testname + ".log"], check=True)

if failed:
    sys.stderr.write("\nAt least one convergence rate deviated from the expected order.\n")
    sys.exit(1)

sys.exit(0)
