#!/usr/bin/env python3

from math import *
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write('Please provide a single argument <testname> to the script\n')
    sys.exit(1)

testname = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]

# remove the old log file
subprocess.call(['rm', testname + '.log'])
print("Removed old log file ({})!".format(testname + '.log'))

# do the runs with different refinement
for i in [0, 1, 2, 3, 4, 5]:
    subprocess.call(['./' + testname] + testargs + ['-Grid.Refinement', str(i),
                                      '-Problem.Name', testname])

# check the rates and append them to the log file
logfile = open(testname + '.log', "r+")

error = []
hmax = []
for line in logfile:
    line = line.strip("\n")
    line = line.strip("\[ConvergenceTest\]")
    line = line.split()
    error.append(float(line[2]))
    hmax.append(float(line[5]))

results = []
logfile.truncate(0)
logfile.write("n\thmax\t\terror\t\trate\n")
logfile.write("-"*50 + "\n")
for i in range(len(error)-1):
    if isnan(error[i]) or isinf(error[i]):
        continue
    if not (error[i] < 1e-12 or error[i+1] < 1e-12):
        rate = (log(error[i])-log(error[i+1]))/(log(hmax[i])-log(hmax[i+1]))
        message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, hmax[i], error[i], rate)
        logfile.write(message)
        results.append(rate)
    else:
        logfile.write("error: exact solution!?")
i = len(error)-1
message = "{}\t{:0.4e}\t{:0.4e}\n\n\n".format(i, hmax[i], error[i])
logfile.write(message)

logfile.close()
print("\nComputed the following convergence rates for {}:\n".format(testname))
subprocess.call(['cat', testname + '.log'])

# check the rates, we expect rates around 2
for r in results:
    if int(round(r)) != 2:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)

sys.exit(0)
