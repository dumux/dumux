#!/usr/bin/env python

from math import *
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt

# this convergence test is for one fixed aperture
if (len(sys.argv) != 3):
    sys.stderr.write("Please provide:\n\
                      1 - the aperture you chose for the simulations\n\
                      2 - 1 or 0 to specify if you want to move the matrix grid points with the aperture or not\n")
    sys.exit(1)

# the chosen aperture
a = sys.argv[1]

# name of the .geo file depending on if the points are to be moved or not
if int(sys.argv[2]) == 1:
    print "\nATTENTION!\nYou specified the grid to be moved with the fracture aperture. Make sure to set MovePoints = true in the input file!\n";
    geoFileName = 'grids/singlefracturequadrilateral_moved.geo'
else:
    geoFileName = 'grids/singlefracturequadrilateral.geo'

# array of permeability ratios and respective markers
k = [1e-4, 1, 1e4]
markers = ['o', '^', 'x']

# loop over the different permeability ratios
for permIndex in range(0, 3):
    # remove the old log file
    subprocess.call(['rm', 'test_1dfracture.log'])
    print("Removed old log file!")

    # for each given number of cells, create new .geo file, mesh & run simulation
    NumCells = np.linspace(10, 100, 10)

    for cells in NumCells:
        tmpGeoFile = open('grids/tmp.geo', "w")
        tmpGeoFile.write("numElements = " + str(int(cells)) + ";\n")
        tmpGeoFile.write("a = " + str(a) + ";\n")

        # copy the rest of the old geo file
        geoFile = open(geoFileName, "r")
        lineCounter = 0
        for line in geoFile:
            # skip the first line (we changed the num elements & aperture)
            if lineCounter > 1:
                tmpGeoFile.write(line)
            lineCounter = lineCounter + 1

        tmpGeoFile.close()
        geoFile.close()

        subprocess.call(['gmsh', '-2', 'grids/tmp.geo'])
        subprocess.call(['./test_analytical', '-Grid.File', 'grids/tmp.msh',
                                              '-Grid.NumCells', str(int(cells)),
                                              '-SpatialParams.FractureAperture', str(a),
                                              '-SpatialParams.FracturePermeability', str(k[permIndex])])

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
    logfile.write("n\ta/L\t\terror_p\t\terror_q\t\trate_p\t\trate_q\n")
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
    plt.figure(1)
    plt.loglog(matrixEps, matrixErrorP, label=r'$k_f = {}$'.format(k[permIndex]), c='b', marker=markers[permIndex])
    plt.loglog(matrixEps, pref, 'b--')
    plt.xlabel(r'$\epsilon_h$', fontsize=22)
    plt.ylabel(r'$eu_{m, L_2}$', fontsize=22)
    plt.legend()

    plt.figure(2)
    plt.loglog(matrixEps, matrixErrorQ, label=r'$k_f = {}$'.format(k[permIndex]), c='r', marker=markers[permIndex])
    plt.loglog(matrixEps, qref, 'r--')
    plt.xlabel(r'$\epsilon_h$', fontsize=22)
    plt.ylabel(r'$eq_{m, L_2}$', fontsize=22)
    plt.legend()

    # quadratically converging reference plots
    p_fref = []
    q_fref = []
    for i in range(0, len(fractureEps)):
        p_fref.append(pow(10, log10(fractureErrorP[0]) - 2*(log10(fractureEps[i]) - log10(fractureEps[0]))))
        q_fref.append(pow(10, log10(fractureErrorQ[0]) - 2*(log10(fractureEps[i]) - log10(fractureEps[0]))))

    # plot fracture results
    plt.figure(3)
    plt.loglog(fractureEps, fractureErrorP, label=r'$k_f = {}$'.format(k[permIndex]), c='b', marker=markers[permIndex])
    plt.loglog(fractureEps, p_fref, 'b--')
    plt.xlabel(r'$\epsilon_h$', fontsize=22)
    plt.ylabel(r'$eu_{f, L_2}$', fontsize=22)
    plt.legend()

    plt.figure(4)
    plt.loglog(fractureEps, fractureErrorQ, label=r'$k_f = {}$'.format(k[permIndex]), c='r', marker=markers[permIndex])
    plt.loglog(fractureEps, q_fref, 'r--')
    plt.xlabel(r'$\epsilon_h$', fontsize=22)
    plt.ylabel(r'$eq_{f, L_2}$', fontsize=22)
    plt.legend()

# show the plots for all permeability ratios
plt.show()
