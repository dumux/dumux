#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


from math import *
import subprocess
import sys
import os
import json

if len(sys.argv) < 3:
    sys.stderr.write("Please provide the following arguments:\n"  \
                     "   - the name of the executable\n" \
                     "   - the name to be used for generated output files\n" \
                     "   - (optional) --expected-rate <rate>, defaults to 2\n" \
                     "   - (optional) runtime arguments to be passed to the executable\n")
    sys.exit(1)

args = [str(i) for i in sys.argv][1:]

# the order of a scheme is a property of the scheme, so the rate to expect is passed in
expectedRate = 2.0
if "--expected-rate" in args:
    i = args.index("--expected-rate")
    expectedRate = float(args[i+1])
    del args[i:i+2]

executableName = args[0]
testName = args[1]
testArgs = args[2:] if len(args) > 2 else ['params.input']

# remove the old log files
if os.path.exists(testName + '.log'):
    subprocess.call(['rm', testName + '.log'])
    print("Removed old log file ({})!".format(testName + '.log'))

# remove the old json files
if os.path.exists(testName + '.json'):
    subprocess.call(['rm', testName + '.json'])
    print("Removed old json file ({})!".format(testName + '.json'))

# do the runs with different refinement
for i in [0, 1, 2, 3]:
    subprocess.call(['./' + executableName]
                    + testArgs
                    + ['-Problem.Name', testName]
                    + ['-Grid.Refinement', str(i)])

def computeRates():
    # check the rates and append them to the log file
    with open(testName + '.json', 'r') as jsonFile:
        data = json.load(jsonFile)
        errorP = data["L2errors"]
        numDofs = data["numDofs"]

    resultsP = []
    with open(testName + '.log', "w") as logFile:
        logFile.truncate(0)
        logFile.write("\nConvergence rates computed for {}:\n\n".format(testName))
        logFile.write("n\terrorP\t\trateP\n")
        logFile.write("-"*50 + "\n")
        for i in range(len(errorP)-1):
            if isnan(errorP[i]) or isinf(errorP[i]):
                continue
            if not ((errorP[i] < 1e-12 or errorP[i+1] < 1e-12)):
                rateP = (log(errorP[i])-log(errorP[i+1]))/log(sqrt(numDofs[i+1]/numDofs[i]))
                message = "{}\t{:0.4e}\t{:0.4e}\n".format(i, errorP[i], rateP)
                logFile.write(message)
                resultsP.append(rateP)
            else:
                logFile.write("error: exact solution!?")
        i = len(errorP)-1
        message = "{}\t{:0.4e}\n".format(i, errorP[i], "")
        logFile.write(message)

    subprocess.call(['cat', testName + '.log'])

    return {"p" : resultsP}

rates = computeRates()

def mean(numbers):
    return float(sum(numbers)) / len(numbers)

# check the rates against the order of the scheme, allowing for the usual preasymptotic slack
if mean(rates["p"]) < 0.9*expectedRate:
    sys.stderr.write("*"*70 + "\n"
                     + "The convergence rates for pressure were not close enough to {}! Test failed.\n".format(expectedRate)
                     + "*"*70 + "\n")
    sys.exit(1)
