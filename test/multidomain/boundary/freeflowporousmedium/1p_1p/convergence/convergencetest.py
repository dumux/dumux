#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later


from math import *
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write('Please provide a single argument <testname> to the script\n')
    sys.exit(1)

executableName = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]
testname = testargs[testargs.index('-Problem.TestCase')+1]
freeflowLogName = testname + "_" + testargs[testargs.index('-FreeFlow.Problem.Name')+1] + '.log'
darcyLogName = testname + "_" + testargs[testargs.index('-Darcy.Problem.Name')+1] + '.log'

if '-Rates' in testargs:
    idx = testargs.index('-Rates')
    rates = testargs[idx+1]
    rates = [float(s) for s in rates.split(',')]
    testargs[idx:idx+2:1] = ""
else:
    rates = [1.7,2.05,1.95,2.05,1.95,2.05]

# remove the old log files
subprocess.call(['rm', freeflowLogName])
print("Removed old log file ({})!".format(freeflowLogName))
subprocess.call(['rm', darcyLogName])
print("Removed old log file ({})!".format(darcyLogName))

# do the runs with different refinement
for i in [0, 1, 2]:
    subprocess.call(['./' + executableName] + testargs + ['-Grid.Refinement', str(i)])

def checkRatesFreeFlow():
    # check the rates and append them to the log file
    with open(freeflowLogName, "r+") as logfile:

        errorP, errorV = [], []
        for line in logfile:
            line = line.strip("\n").strip("\[ConvergenceTest\]").split()
            errorP.append(float(line[2]))
            errorV.append(float(line[5]))

        resultsP, resultsV = [], []
        logfile.truncate(0)
        logfile.write("n\terrorP\t\trateP\t\terrorV\t\trateV\n")
        logfile.write("-"*50 + "\n")
        for i in range(len(errorP)-1):
            if isnan(errorP[i]) or isinf(errorP[i]):
                continue
            if not ((errorP[i] < 1e-12 or errorP[i+1] < 1e-12) and (errorV[i] < 1e-12 or errorV[i+1] < 1e-12)):
                rateP = (log(errorP[i])-log(errorP[i+1]))/log(2)
                rateV = (log(errorV[i])-log(errorV[i+1]))/log(2)
                message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, errorP[i], rateP, errorV[i], rateV)
                logfile.write(message)
                resultsP.append(rateP)
                resultsV.append(rateV)
            else:
                logfile.write("error: exact solution!?")
        i = len(errorP)-1
        message = "{}\t{:0.4e}\t\t{}\t{:0.4e}\t\t{}\n".format(i, errorP[i], "", errorV[i], "")
        logfile.write(message)

    print("\nComputed the following convergence rates for {}:\n".format(testname))

    subprocess.call(['cat', freeflowLogName])

    return {"p" : resultsP, "v" : resultsV}

def checkRatesDarcy():
    # check the rates and append them to the log file
    with open(darcyLogName, "r+") as logfile:

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

    print("\nComputed the following convergence rates for {}:\n".format(testname))

    subprocess.call(['cat', darcyLogName])

    return {"p" : resultsP}

def checkRatesFreeFlowAndDarcy():
    resultsFreeFlow = checkRatesFreeFlow()
    resultsDarcy = checkRatesDarcy()

    def mean(numbers):
        return float(sum(numbers)) / len(numbers)

    # check the rates, we expect rates around 2
    if mean(resultsFreeFlow["p"]) < rates[1] and mean(resultsFreeFlow["p"]) < rates[0]:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for pressure were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)

    if mean(resultsFreeFlow["v"]) < rates[3] and mean(resultsFreeFlow["v"]) < rates[2]:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for velocity were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)

    if mean(resultsDarcy["p"]) < rates[5] and mean(resultsDarcy["p"]) < rates[4]:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for pressure were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)


checkRatesFreeFlowAndDarcy()
