#!/usr/bin/env python3

from math import *
import subprocess
import sys
import os

if len(sys.argv) < 3:
    sys.stderr.write("Please provide the following arguments:\n"  \
                     "   - the name of the executable\n" \
                     "   - the name to be used for generated output files\n" \
                     "   - (optional) runtime arguments to be passed to the executable\n")
    sys.exit(1)

executableName = sys.argv[1]
testName = sys.argv[2]
testArgs = [str(i) for i in sys.argv][3:] if len(sys.argv) > 3 else ['params.input']

# remove the old log files
if os.path.exists(testName + '.log'):
    subprocess.call(['rm', testName + '.log'])
    print("Removed old log file ({})!".format(testName + '.log'))

# do the runs with different refinement
for i in [0, 1, 2, 3]:
    subprocess.call(['./' + executableName]
                    + testArgs
                    + ['-Problem.Name', testName]
                    + ['-Grid.Refinement', str(i)])

def computeRates():
    # check the rates and append them to the log file
    logFile = open(testName + '.log', "r+")

    errorP = []
    for line in logFile:
        line = line.strip("\n")
        line = line.strip("\[ConvergenceTest\]")
        line = line.split()
        errorP.append(float(line[2]))

    resultsP = []
    logFile.truncate(0)
    logFile.write("\nConvergence rates computed for {}:\n\n".format(testName))
    logFile.write("n\terrorP\t\trateP\n")
    logFile.write("-"*50 + "\n")
    for i in range(len(errorP)-1):
        if isnan(errorP[i]) or isinf(errorP[i]):
            continue
        if not ((errorP[i] < 1e-12 or errorP[i+1] < 1e-12)):
            rateP = (log(errorP[i])-log(errorP[i+1]))/log(2)
            message = "{}\t{:0.4e}\t{:0.4e}\n".format(i, errorP[i], rateP)
            logFile.write(message)
            resultsP.append(rateP)
        else:
            logFile.write("error: exact solution!?")
    i = len(errorP)-1
    message = "{}\t{:0.4e}\n".format(i, errorP[i], "")
    logFile.write(message)
    logFile.close()

    subprocess.call(['cat', testName + '.log'])

    return {"p" : resultsP}

rates = computeRates()

def mean(numbers):
    return float(sum(numbers)) / len(numbers)

# check the rates, we expect rates around 2
if mean(rates["p"]) < 1.8:
    sys.stderr.write("*"*70 + "\n" + "The convergence rates for pressure were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
    sys.exit(1)
