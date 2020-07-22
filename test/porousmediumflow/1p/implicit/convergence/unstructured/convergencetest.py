#!/usr/bin/env python3

from math import *
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write('Please provide a single argument <testname> to the script\n')
    sys.exit(1)

executableName = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]
testname = testargs[2]

# remove the old log files
subprocess.call(['rm', testname + '.log'])
print("Removed old log file ({})!".format(testname + '.log'))

# do the runs with different refinement
for i in [0, 1, 2, 3]:
    subprocess.call(['./' + executableName] + testargs + ['-Grid.Refinement', str(i)])

def checkRates():
    # check the rates and append them to the log file
    logfile = open(testname + '.log', "r+")

    errorP = []
    for line in logfile:
        line = line.strip("\n")
        line = line.strip("\[ConvergenceTest\]")
        line = line.split()
        errorP.append(float(line[2]))

    resultsP = []
    logfile.truncate(0)
    logfile.write("n\terrorP\t\trateP\n")
    logfile.write("-"*50 + "\n")
    for i in range(len(errorP)-1):
        if isnan(errorP[i]) or isinf(errorP[i]):
            continue
        if not ((errorP[i] < 1e-12 or errorP[i+1] < 1e-12)):
            rateP = (log(errorP[i])-log(errorP[i+1]))/log(2)
            message = "{}\t{:0.4e}\t{:0.4e}\n".format(i, errorP[i], rateP)
            logfile.write(message)
            resultsP.append(rateP)
        else:
            logfile.write("error: exact solution!?")
    i = len(errorP)-1
    message = "{}\t{:0.4e}\n".format(i, errorP[i], "")
    logfile.write(message)

    logfile.close()
    print("\nComputed the following convergence rates for {}:\n".format(testname))

    subprocess.call(['cat', testname + '_darcy.log'])

    return {"p" : resultsP}

def checkConvRates():
    rates = checkRates()

    def mean(numbers):
        return float(sum(numbers)) / len(numbers)

    # check the rates, we expect rates around 2
    if mean(rates["p"]) < 2.05 and mean(rates["p"]) < 1.95:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for pressure were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)


checkConvRates()
