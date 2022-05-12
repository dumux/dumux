#!/usr/bin/env python3

import os
import sys
import subprocess
from math import pow
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt

def initializePlot(pltIdx=1):
    plt.figure(pltIdx)
    plt.clf()

def plotResultsFromSubFolder(fileName, idx, label, pltIdx=1):
    plt.figure(pltIdx)
    errors = np.genfromtxt(fileName, delimiter=',',skip_header=0, dtype=float)
    h = errors[:,idx[0]]
    values = errors[:,idx[1]]
    plt.loglog(h, values, label=label, marker='o')
    plt.legend()
    plt.grid(True)


def plotReferenceCurve(fileName, idx, pltIdx=1):
    plt.figure(pltIdx)
    errors = np.genfromtxt(fileName, delimiter=',',skip_header=0, dtype=float)
    h = errors[:,idx[0]]
    values = errors[:,idx[1]]
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


def runConvergenceTest(exeName, outputName, paramsFile, numRefinements=5, extraArgs=[]):
    for i in range(numRefinements):
        baseCall = [
            "./" + exeName, paramsFile,
            "-Grid.Refinement", str(i),
            "-Vtk.OutputName", str(outputName)
        ] + extraArgs
        subprocess.run(" ".join(str(x) for x in baseCall), shell=True, executable='/bin/bash')

def runConvergenceTestAndMoveToFolder(
    exeName, folderName, outputName, paramsFile, numRefinements=5, extraArgs=[]
):
    os.makedirs(folderName, exist_ok=True)
    runConvergenceTest(exeName, outputName, paramsFile, numRefinements, extraArgs)
    moveResultsToSubFolder(folderName)

def runTests(
    testNames, paramsFile, testRuns, subRuns, numRefinements
):
    removePreviousResults()
    for l in range(len(testNames)):
        testName = testNames[l]
        exes = testRuns[l]
        for k in range(len(exes)):
            exe = exes[k]
            compile(exe[0])
            testargs = exe[1]
            folderName = testName
            subRun = subRuns[l][k]
            for j in range(len(subRun)):
                run = subRun[j]
                outputName = run[0][0] + "-" + run[0][1]
                extraArgs = testargs + run[1] + ["-Problem.Name " + outputName]
                runConvergenceTestAndMoveToFolder(
                    exe[0], folderName, outputName, paramsFile, numRefinements, extraArgs
                )

def plotResults(
    testNames, testRuns, subRuns, indices
):
    for l in range(len(testNames)):
        for m in range(len(indices)):
            initializePlot(m+1)
        testName = testNames[l]
        exes = testRuns[l]
        for k in range(len(exes)):
            folderName = testName
            subRun = subRuns[l][k]
            for j in range(len(subRun)):
                run = subRun[j]
                outputName = run[0][0] + "-" + run[0][1]
                fileName = folderName + "/" + outputName + ".csv"
                for m in range(len(indices)):
                    plotResultsFromSubFolder(fileName, indices[m][0], run[0][1]+ "_" + indices[m][1],m+1)

        for m in range(len(indices)):
            plotReferenceCurve(fileName, indices[m][0], m+1)
            setAxesLabels("$h_\mathrm{ref}$", "$e_{" + indices[m][1] + "}$")
            savePlot(folderName + "/" + indices[m][1] + "-error-" + folderName + ".pdf")
