#!/usr/bin/env python3

import os
import sys
import subprocess
from math import pow
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt

def initializePlot():
    plt.figure()
    plt.clf()


def plotResultsFromSubFolder(fileName, idx, label):
    errors = np.genfromtxt(fileName, delimiter=',',skip_header=0, dtype=float)
    h = errors[:,2]
    values = errors[:,idx]
    plt.loglog(h, values, label=label, marker='o')
    plt.legend()
    plt.grid(True)


def plotReferenceCurve(fileName, idx):
    errors = np.genfromtxt(fileName, delimiter=',',skip_header=0, dtype=float)
    h = errors[:,2]
    values = errors[:,idx]
    referenceValue = values[0]/1.2
    plt.loglog(h,
              [referenceValue*h[i]*h[i] / (h[0]*h[0]) for i in range(len(values))],
              marker="*",
              linestyle="--",
              label="$\mathcal{O}(h^2)$")
    plt.legend()


def setAxesLabels(xlabel, ylabel):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


def savePlot(name):
    plt.savefig(name, bbox_inches="tight")

def exitWithError(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def fileHasExtension(file, ext):
    fileExt = os.path.splitext(file)[1]
    return fileExt == ext or fileExt.strip(".") == ext


def findFilesWithExtension(ext):
    return list(filter(
        lambda f: os.path.isfile(f) and fileHasExtension(f, ext),
        os.listdir()
    ))


def applyToFilesWithExtension(ext, functor):
    for file in findFilesWithExtension(ext):
        functor(file)


def removeFilesWithExtension(ext):
    applyToFilesWithExtension(ext, os.remove)


def moveFilesWithExtension(ext, targetFolder):
    def makeTargetPath(f):
        return os.path.join(targetFolder, f)
    applyToFilesWithExtension(ext, lambda f: os.replace(f, makeTargetPath(f)))


def removePreviousResults():
    removeFilesWithExtension("vtu")
    removeFilesWithExtension("vtp")
    removeFilesWithExtension("pvd")
    removeFilesWithExtension("log")
    removeFilesWithExtension("csv")


def moveResultsToSubFolder(folderName):
    moveFilesWithExtension("vtu", folderName)
    moveFilesWithExtension("vtp", folderName)
    moveFilesWithExtension("pvd", folderName)
    moveFilesWithExtension("log", folderName)
    moveFilesWithExtension("csv", folderName)

def compile(exeName):
    subprocess.run(["make", exeName], check=True)


def runConvergenceTest(exeName, outputName, numRefinements=5, extraArgs=[]):
    for i in range(numRefinements):
        baseCall = [
            "./" + exeName, "params.input",
            "-Grid.Refinement", str(i),
            "-Vtk.OutputName", str(outputName)
        ] + extraArgs
        subprocess.run(" ".join(str(x) for x in baseCall), shell=True, executable='/bin/bash')


def runConvergenceTestAndMoveToFolder(
    exeName, folderName, outputName, numRefinements=5, extraArgs=[]
):
    os.makedirs(folderName, exist_ok=True)
    runConvergenceTest(exeName, outputName, numRefinements, extraArgs)
    moveResultsToSubFolder(folderName)

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Run script for the convergence test."
    )
    parser.add_argument("-n", "--num-refinements", required=False, default=5)
    #parser.add_argument("-e", "--executable", type=str, required=True)
    #parser.add_argument("-o", "--outputname", type=str, required=True)
    args = vars(parser.parse_args())

    numRefinements = int(args["num_refinements"])
    #exe = str(args["executable"])
    #outputName = str(args["outputname"])

    removePreviousResults()

    testPrefix = "donea-momentum"
    executables = ["test_ff_stokes_donea_momentum_diamond_simplex"]
    gridArgs = [ [], ["-Grid.File ./grids/unstructured_simplex.msh"]]
    gridNames = ["simplex-structured", "simplex-unstructured" ]

    ##################################################
    # structured simplex grid
    ##################################################

    for exe in executables:
        compile(exe)
        for i in range(len(gridArgs)):
            initializePlot()
            # weak symmetry
            folderName = testPrefix + "-" + gridNames[i]

            methods = ["weak-sym", "unsym"]
            methodArgs = [ [], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]]

            for j in range(len(methodArgs)):
                method = methods[j]
                outputName = folderName + "-" + method
                extraArgs = gridArgs[i] + methodArgs[j]
                runConvergenceTestAndMoveToFolder(
                    exe, folderName, outputName, numRefinements, extraArgs
                )
                fileName = folderName + "/" + "errors.csv"
                plotResultsFromSubFolder(fileName, 3, method)

            plotReferenceCurve(fileName, 3)
            setAxesLabels("$h_\mathrm{max}$", "$e_{v}$")
            savePlot(folderName + "/" + "vel-error-" + folderName + ".pdf")