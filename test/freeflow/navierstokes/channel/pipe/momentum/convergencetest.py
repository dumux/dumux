#!/usr/bin/env python3
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


from math import *
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write("Please provide a single argument <testname> to the script\n")
    sys.exit(1)

testname = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]

# remove the old log file
subprocess.call(["rm", testname + ".log"])
print("Removed old log file ({})!".format(testname + ".log"))

# do the runs with different refinement
for i in [0, 1, 2]:
    subprocess.call(
        ["./" + testname] + testargs + ["-Grid.Refinement", str(i), "-Problem.Name", testname]
    )

# check the rates and append them to the log file
logfile = open(testname + ".log", "r+")

errorVx = []
errorVy = []
for line in logfile:
    line = line.strip("\n")
    line = line.strip("\[ConvergenceTest\]")
    line = line.split()
    print(line)
    errorVx.append(float(line[5]))
    errorVy.append(float(line[8]))

resultsVx = []
resultsVy = []
logfile.truncate(0)
logfile.write("n\terrorVx\t\trateVx\t\terrorVy\t\trateVy\n")
logfile.write("-" * 50 + "\n")
for i in range(len(errorVx) - 1):
    if isnan(errorVx[i]) or isinf(errorVx[i]):
        continue
    if not (
        (errorVx[i] < 1e-12 or errorVx[i + 1] < 1e-12)
        and (errorVy[i] < 1e-12 or errorVy[i + 1] < 1e-12)
    ):
        rateVx = (log(errorVx[i]) - log(errorVx[i + 1])) / log(2)
        rateVy = (log(errorVy[i]) - log(errorVy[i + 1])) / log(2)
        message = f"{i}\t{errorVx[i]:0.4e}\t{rateVx:0.4e}\t{errorVy[i]:0.4e}\t{rateVy:0.4e}\n"
        logfile.write(message)
        resultsVx.append(rateVx)
        resultsVy.append(rateVy)
    else:
        logfile.write("error: exact solution!?")
i = len(errorVx) - 1
message = f"{i}\t{errorVx[i]:0.4e}\t\t\t{errorVy[i]:0.4e}\t\t\n"
logfile.write(message)

logfile.close()
print("\nComputed the following convergence rates for {}:\n".format(testname))
subprocess.call(["cat", testname + ".log"])


def mean(numbers):
    return float(sum(numbers)) / len(numbers)


if mean(resultsVx) < 2.05 and mean(resultsVx) < 1.9 and max(errorVx) > 1e-15:
    sys.stderr.write(
        "*" * 70
        + "\n"
        + "The convergence rates for x-velocity were not close enough to 2! Test failed.\n"
        + "*" * 70
        + "\n"
    )
    sys.exit(1)

if mean(resultsVy) < 2.05 and mean(resultsVy) < 1.9 and max(errorVy) > 1e-15:
    sys.stderr.write(
        "*" * 70
        + "\n"
        + "The convergence rates for y-velocity were not close enough to 2! Test failed.\n"
        + "*" * 70
        + "\n"
    )
    sys.exit(1)

sys.exit(0)
