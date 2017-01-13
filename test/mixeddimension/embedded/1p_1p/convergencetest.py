#!/usr/bin/env python

from math import *
import subprocess
import sys

if len(sys.argv) != 2:
    sys.err.write('Please provide a single argument <testname> to the script')

testname = str(sys.argv[1])

subprocess.call(['rm', testname + '.log'])
print("Removed old log file ({})!".format(testname + '.log'))
for i in [10, 20, 30, 40, 50, 60]:
    subprocess.call(['./' + testname, '-ParameterFile', 'test_1p1dstokes.input',
                                      '-TissueGrid.Cells', 3*(str(i)+" "),
                                      '-VesselGrid.Cells', str(i),
                                      '-Problem.Name', testname])

logfile = open(testname + '.log', "r+")

bulkError = []
bulkHmax = []
lowDimError = []
lowDimHmax = []
for line in logfile:
    line = line.strip("\n")
    line = line.split()
    bulkHmax.append(float(line[0]))
    lowDimHmax.append(float(line[1]))
    bulkError.append(float(line[2]))
    lowDimError.append(float(line[3]))

logfile.truncate(0)
for name, hmax, error in (("3D", bulkHmax, bulkError), ("1D", lowDimHmax, lowDimError)):
    logfile.write("[" + name + "]\thmax\t\terror\t\trate\n")
    logfile.write("-"*50 + "\n")
    for i in range(len(error)-1):
        if isnan(error[i]) or isinf(error[i]):
            continue
        if not (error[i] < 1e-12 or error[i] < 1e-12):
            rate = (log(error[i])-log(error[i+1]))/(log(hmax[i])-log(hmax[i+1]))
            message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, hmax[i], error[i], rate)
            logfile.write(message)
        else:
            logfile.write("error: exact solution!?")
    i = len(error)-1
    message = "{}\t{:0.4e}\t{:0.4e}\n\n\n".format(i, hmax[i], error[i])
    logfile.write(message)

logfile.close()
print("\nComputed the following convergence rates:\n")
subprocess.call(['cat', testname + '.log'])
