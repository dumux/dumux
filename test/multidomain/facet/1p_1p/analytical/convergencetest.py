#!/usr/bin/env python3

from math import *
import os.path
import subprocess
import sys

if len(sys.argv) < 3:
    sys.stderr.write("Invalid number of arguments given. Please provide the following arguments (in this order):\n\
                      - the name of the executable (either test_facetcoupling_tpfa_1p1p or test_facetcoupling_box_1p1p)\n\
                      - the fracture permeabilities you want to perform the convergence test for")
    sys.exit(1)

# name of the executable
execName = sys.argv[1]

# array of permeabilities
k = []
for i in range(2, len(sys.argv)):
    k.append(sys.argv[i])

# loop over the different permeability ratios
for permIndex in range(0, len(k)):
    # check whether executable exists
    if not os.path.isfile(execName):
        sys.exit(77)

    # remove the old log file
    subprocess.call(['rm', execName + '.log'])
    print("Removed old log file (" + execName + '.log' + ")!")

    # for each given number of cells, create new .geo file, mesh & run simulation
    numCells = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for cells in numCells:
        geoFile = open("grids/hybridgrid.geo", 'r')
        lines = geoFile.readlines()
        firstLine = lines[0]
        firstLine = firstLine.split(" ")
        if firstLine[0] != "numElemsPerSide":
            sys.stderr.write("Geo file is expected to have - numElemsPerSide = xx - in the first line!\n")
            sys.exit(1)

        # write a temporary geo-file using num elements
        tmpGeoFile = open("grids/" + execName + ".geo", 'w')
        lines[0] = "numElemsPerSide = " + str(int(cells)) + ";\n"
        for line in lines:
            tmpGeoFile.write(line)

        geoFile.close()
        tmpGeoFile.close()

        subprocess.call(['gmsh', '-format', 'msh2', '-2', "grids/" + execName + ".geo"])
        subprocess.call(['./' + execName, 'params.input',
                                          '-Vtk.OutputName', execName,
                                          '-Grid.File', "grids/" + execName + ".msh",
                                          '-Grid.NumElemsPerSide', str(int(cells)),
                                          '-LowDim.SpatialParams.Permeability', str(k[permIndex]),
                                          '-Problem.OutputFileName', execName + '.log'])

        # remove geo and msh file
        subprocess.call(['rm', "grids/" + execName + ".geo"])
        subprocess.call(['rm', "grids/" + execName + ".msh"])

    # check the rates and append them to the log file
    logfile = open(execName + '.log', "r+")

    eps = []
    l2_matrix = []
    l2_fracture = []

    for line in logfile:
        line = line.strip("\n")
        line = line.split(",")
        eps.append(float(line[0]))
        l2_matrix.append(float(line[1]))
        l2_fracture.append(float(line[3]))

    matrix_rates = []
    fracture_rates = []

    logfile.truncate(0)
    logfile.write("Matrix domain:\n")
    logfile.write("-"*50 + "\n")
    logfile.write("n\ta/L\t\terror_p\t\trate_p\n")
    logfile.write("-"*50 + "\n")
    for i in range(len(l2_matrix)-1):
        if isnan(l2_matrix[i]) or isinf(l2_matrix[i]):
            sys.stderr.write("l2 error norm is not a number!\n")
            sys.exit(1)
        if not (l2_matrix[i] < 1e-12 or l2_matrix[i+1] < 1e-12):
            rateP = (log(l2_matrix[i])-log(l2_matrix[i+1]))/(log(eps[i])-log(eps[i+1]))
            message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, eps[i], l2_matrix[i], rateP)
            logfile.write(message)
            matrix_rates.append(rateP)
        else:
            sys.stderr.write("error < 1e-12 does not seem reasonable!\n")
            sys.exit(1)

    # test is failed if the rate in the matrix is lower than 1.9
    if abs(matrix_rates[len(matrix_rates)-1]) < 1.9:
        sys.stderr.write("\n\nComputed a convergence rate of " + str(abs(matrix_rates[len(matrix_rates)-1])) + ", which is below the treshold of 1.9. Test fails...\n")
        sys.exit(1)

    # write last output for the matrix
    i = len(l2_matrix)-1
    message = "{}\t{:0.4e}\t{:0.4e}\n".format(i, eps[i], l2_matrix[i])
    logfile.write(message)

    logfile.write("\n\nFracture domain:\n")
    logfile.write("-"*50 + "\n")
    logfile.write("n\ta/L\t\terror_p\t\trate_p\n")
    logfile.write("-"*50 + "\n")
    for i in range(len(l2_fracture)-1):
        if isnan(l2_fracture[i]) or isinf(l2_fracture[i]):
            sys.stderr.write("l2 error norm is not a number!\n")
            sys.exit(1)
        if not (l2_fracture[i] < 1e-12 or l2_fracture[i] < 1e-12):
            rateP = (log(l2_fracture[i])-log(l2_fracture[i+1]))/(log(eps[i])-log(eps[i+1]))
            message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, eps[i], l2_fracture[i], rateP)
            logfile.write(message)
            fracture_rates.append(rateP)
        else:
            sys.stderr.write("error < 1e-12 does not seem reasonable!\n")
            sys.exit(1)

    # write last output for the matrix
    i = len(l2_fracture)-1
    message = "{}\t{:0.4e}\t{:0.4e}\n".format(i, eps[i], l2_fracture[i])
    logfile.write(message)
    logfile.close()

    print("\nComputed the following convergence rates for:")
    subprocess.call(['cat', execName + '.log'])
