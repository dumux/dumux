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

def runTestsAndPlotResults(
    testPrefix, executables, grids, methods, numRefinements
):
    removePreviousResults()

    for k in range(len(executables)):
        exe = executables[k]
        compile(exe)
        gridsExe = grids[k]
        for i in range(len(gridsExe)):
            grid = gridsExe[i]
            initializePlot()
            # weak symmetry
            folderName = testPrefix + "-" + grid[0]

            for j in range(len(methods)):
                method = methods[j]
                outputName = folderName + "-" + method[0]
                extraArgs = grid[1] + method[1]
                runConvergenceTestAndMoveToFolder(
                    exe, folderName, outputName, numRefinements, extraArgs
                )
                fileName = folderName + "/" + "errors.csv"
                plotResultsFromSubFolder(fileName, 3, method[0])

            plotReferenceCurve(fileName, 3)
            setAxesLabels("$h_\mathrm{ref}$", "$e_{v}$")
            savePlot(folderName + "/" + "vel-error-" + folderName + ".pdf")
