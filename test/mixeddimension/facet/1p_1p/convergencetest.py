#!/usr/bin/env python

from math import *
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# remove the old log file
subprocess.call(['rm', 'test_1dfracture.log'])
print("Removed old log file!")

# do the runs with different grids
# this convergence test is for a fixed aperture
a = 1e-5;

# for each given number of cells, create new .geo file, mesh & run simulation
NumCells = []
for i in range(0,19):
    NumCells.append(10 + i*5)

for cells in NumCells:
    tmpGeoFile = open('grids/tmp.geo', "w")
    tmpGeoFile.write("numElements = " + str(int(cells)) + ";\n")

    # copy the rest of the old geo file
    geoFile = open('grids/singlefracturequadrilateral.geo', "r")
    lineCounter = 0
    for line in geoFile:
        # skip the first line (we changed the num elements)
        if lineCounter != 0:
            tmpGeoFile.write(line)
        lineCounter = lineCounter + 1

    tmpGeoFile.close()
    geoFile.close()

    subprocess.call(['gmsh', '-2', 'grids/tmp.geo'])
    subprocess.call(['./test_analytical', '-Grid.File', 'grids/tmp.msh',
                                          '-Grid.NumCells', str(cells),
                                          '-SpatialParams.FractureAperture', str(a)])

    # remove geo and msh file
    subprocess.call(['rm', 'grids/tmp.geo'])
    subprocess.call(['rm', 'grids/tmp.msh'])

# check the rates and append them to the log file
logfile = open('test_1dfracture.log', "r+")

matrixErrorP = []
matrixErrorQ = []
matrixEps = []

fractureErrorP = []
fractureErrorQ = []
fractureEps = []

for line in logfile:
    line = line.strip("\n")
    line = line.split()
    matrixEps.append(float(line[0]))
    matrixErrorP.append(float(line[1]))
    matrixErrorQ.append(float(line[2]))
    fractureEps.append(float(line[3]))
    fractureErrorP.append(float(line[4]))
    fractureErrorQ.append(float(line[5]))

matrixRatesP = []
matrixRatesQ = []
fractureRatesP = []
fractureRatesQ = []

logfile.truncate(0)
logfile.write("Matrix domain:\n")
logfile.write("-"*50 + "\n")
logfile.write("n\ta/L\t\terror_p\t\terror_1\t\trate_p\t\trate_q\n")
logfile.write("-"*50 + "\n")
for i in range(len(matrixErrorP)-1):
    if isnan(matrixErrorP[i]) or isinf(matrixErrorP[i]) or isnan(matrixErrorQ[i]) or isinf(matrixErrorQ[i]):
        continue
    if not (matrixErrorP[i] < 1e-12 or matrixErrorQ[i] < 1e-12):
        rateP = (log(matrixErrorP[i])-log(matrixErrorP[i+1]))/(log(matrixEps[i])-log(matrixEps[i+1]))
        rateQ = (log(matrixErrorQ[i])-log(matrixErrorQ[i+1]))/(log(matrixEps[i])-log(matrixEps[i+1]))
        message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, matrixEps[i], matrixErrorP[i], matrixErrorQ[i], rateP, rateQ)
        logfile.write(message)
        matrixRatesP.append(rateP)
        matrixRatesQ.append(rateQ)
    else:
        logfile.write("error < 1e-12 seems impossible. Did you use the exact solution!?")

# write last output for the matrix
i = len(matrixErrorP)-1
message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, matrixEps[i], matrixErrorP[i], matrixErrorQ[i])
logfile.write(message)

logfile.write("\n\nFracture domain:\n")
logfile.write("-"*50 + "\n")
logfile.write("n\ta/L\t\terror_p\t\terror_1\t\trate_p\t\trate_q\n")
logfile.write("-"*50 + "\n")
for i in range(len(fractureErrorP)-1):
    if isnan(fractureErrorP[i]) or isinf(fractureErrorP[i]) or isnan(fractureErrorQ[i]) or isinf(fractureErrorQ[i]):
        continue
    if not (fractureErrorP[i] < 1e-12 or fractureErrorQ[i] < 1e-12):
        rateP = (log(fractureErrorP[i])-log(fractureErrorP[i+1]))/(log(fractureEps[i])-log(fractureEps[i+1]))
        rateQ = (log(fractureErrorQ[i])-log(fractureErrorQ[i+1]))/(log(fractureEps[i])-log(fractureEps[i+1]))
        message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, fractureEps[i], fractureErrorP[i], fractureErrorQ[i], rateP, rateQ)
        logfile.write(message)
        fractureRatesP.append(rateP)
        fractureRatesQ.append(rateQ)
    else:
        logfile.write("error < 1e-12 seems impossible. Did you use the exact solution!?")

# write last output for the matrix
i = len(fractureErrorP)-1
message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, fractureEps[i], fractureErrorP[i], fractureErrorQ[i])
logfile.write(message)
logfile.close()

print("\nComputed the following convergence rates for:")
subprocess.call(['cat', 'test_1dfracture.log'])

print("\nplotting...")
# quadratically converging reference plots
pref = []
qref = []
for i in range(0, len(matrixEps)):
    pref.append(pow(10, log10(matrixErrorP[0]) - 2*(log10(matrixEps[i]) - log10(matrixEps[0]))))
    qref.append(pow(10, log10(matrixErrorQ[0]) - 2*(log10(matrixEps[i]) - log10(matrixEps[0]))))

# plot matrix results
plt.loglog(matrixEps, matrixErrorP, label=r'$ep_{L_2}$')
plt.loglog(matrixEps, matrixErrorQ, label=r'$eq_{L_2}$')
plt.loglog(matrixEps, pref, 'r--')
plt.loglog(matrixEps, qref, 'b--')
plt.xlabel(r'$\epsilon_h$', fontsize=18)
plt.ylabel(r'$e_{L_2}$', fontsize=18)
plt.legend()
plt.show()

# plot fracture results
# plt.loglog(fractureEps, fractureErrorP, label=r'$ep_{L_2}$')
# plt.loglog(fractureEps, fractureErrorQ, label=r'$eq_{L_2}$')
# plt.xlabel(r'$\epsilon_h$', fontsize=18)
# plt.ylabel(r'$e_{L_2}$', fontsize=18)
# plt.legend()
# plt.show()
