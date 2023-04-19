# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Helper script to compute discrete l2 error from output files
Note: better do this in the code for increased precision
"""

import argparse
import csv
import sys


def parseCommandLine():
    """Auxiliary function that parses command line arguments"""

    parser = argparse.ArgumentParser(
        prog="python " + sys.argv[0],
        description="Calculate the l2 error of csv data files.",
    )
    parser.add_argument("-f1", "--reference", type=str, required=True, help="Reference csv-file")
    parser.add_argument(
        "-f2", "--newSolution", type=str, required=True, help="NewSolution csv-file"
    )
    parser.add_argument(
        "-xMin",
        "--xMin",
        type=float,
        required=False,
        default=-1e99,
        help="Restrict data to x>xMin",
    )
    parser.add_argument(
        "-xMax",
        "--xMax",
        type=float,
        required=False,
        default=1e99,
        help="Restrict data to x>xMax",
    )
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument("-x1", "--xData1", type=int, help="Column index of x data in reference")
    group1.add_argument("-x1Name", "--xDataName1", type=str, help="Name x data in reference")
    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument("-x2", "--xData2", type=int, help="Column index of x data in newSolution")
    group2.add_argument("-x2Name", "--xDataName2", type=str, help="Name x data in newSolution")
    group3 = parser.add_mutually_exclusive_group(required=True)
    group3.add_argument("-y1", "--yData1", type=int, help="Column index of y data in reference")
    group3.add_argument("-y1Name", "--yDataName1", type=str, help="Name y data in reference")
    group4 = parser.add_mutually_exclusive_group(required=True)
    group4.add_argument("-y2", "--yData2", type=int, help="Column index of y data in newSolution")
    group4.add_argument("-y2Name", "--yDataName2", type=str, help="Name y data in newSolution")
    parser.add_argument("-p", "--percent", action="store_true", help="Print errors in percent")
    parser.add_argument("-f", "--force", action="store_true", help="Ignore 'not-matching' errors")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbosity of the script")
    return vars(parser.parse_args())


def computeIndices(reference, args):
    """Compute the indices of the reference data"""
    if args["xDataName1"] is not None:
        iRefX = reference[0].index(args["xDataName1"])
    else:
        iRefX = args["xData1"]
    if args["yDataName1"] is not None:
        iRefY = reference[0].index(args["yDataName1"])
    else:
        iRefY = args["yData1"]

    if args["xDataName2"] is not None:
        iSolX = reference[0].index(args["xDataName2"])
    else:
        iSolX = args["xData2"]
    if args["yDataName2"] is not None:
        iSolY = reference[0].index(args["yDataName2"])
    else:
        iSolY = args["yData2"]

    return (
        (iRefX, iRefY),
        (iSolX, iSolY),
    )


def readL2ErrorData(args):
    """Compute L2 error: main driver"""

    with open(args["reference"], "rb") as referenceFile:
        reference = list(csv.reader(referenceFile))
    with open(args["newSolution"], "rb") as newSolutionFile:
        newSolution = list(csv.reader(newSolutionFile))

    iRef, iSol = computeIndices(reference, args)

    if reference[0][iRef[0]] != reference[0][iSol[0]] and not args["force"]:
        print(
            "X-Identifier not equal: ref=",
            reference[0][iRef[0]],
            ",new=",
            reference[0][iSol[0]],
            ". Aborting! (Use -f to continue anyway)",
        )
        sys.exit(1)

    if reference[0][iRef[1]] != newSolution[0][iSol[0]] and not args["force"]:
        print(
            "Y-Identifier not equal. ref=",
            reference[0][iRef[1]],
            ",new=",
            newSolution[0][iSol[1]],
            ". Aborting! (Use -f to continue anyway)",
        )
        sys.exit(2)

    if len(reference) != len(newSolution):
        print(
            "Length of reference and newSolution not equal: ref=",
            len(reference),
            ",new=",
            len(newSolution),
            ". Aborting!",
        )
        sys.exit(3)

    return {
        "reference": reference,
        "iRef": iRef,
        "newSolution": newSolution,
        "iSol": iSol,
    }


def computeL2ErrorSquared(args, reference, iRef, newSolution, iSol):
    """Compute L2 error"""
    sumError = 0.0
    sumReference = 0.0
    sumDistance = 0.0
    numPoints = 0

    for i in range(1, len(reference)):
        coordReference = float(reference[i][iRef[0]])
        coordNewSolution = float(newSolution[i][iSol[0]])
        if coordReference != coordNewSolution:
            print(
                "Coordinates not equal: ref=",
                coordReference,
                ",new=",
                coordNewSolution,
                ". Aborting!",
            )
            sys.exit(4)
        if coordReference < float(args["xMin"]) or coordReference > float(args["xMax"]):
            continue

        if i == 1:
            distance = 0.5 * (float(reference[2][iRef[0]]) - float(reference[1][iRef[0]]))
        elif i == len(reference) - 1:
            distA = float(reference[len(reference) - 1][iRef[0]])
            distB = float(reference[len(reference) - 2][iRef[0]])
            distance = 0.5 * (distA - distB)
        else:
            distance = 0.5 * (float(reference[i + 1][iRef[0]]) - float(reference[i - 1][iRef[0]]))
        sumError += (
            (float(reference[i][iRef[1]]) - float(newSolution[i][iSol[1]])) ** 2
        ) * distance
        sumReference += ((float(reference[i][iRef[1]])) ** 2) * distance
        sumDistance += distance
        numPoints += 1

    if numPoints < 999 and not args["force"]:
        print(
            "Warning: numPoints=",
            numPoints,
            " is low, could result in bad the error approximation."
            " (Use -f to suppress this warning)",
        )

    return {"absolute": sumError / sumDistance, "relative": sumError / sumReference}


def printL2Error(args, absolute, relative):
    """Print L2 error"""

    # numPoints is needed, resulting from the equidistant integration
    l2normAbs = (absolute) ** 0.5
    # numPoints cancels out for equidistant integration
    l2normRel = (relative) ** 0.5

    if args["percent"]:
        print(
            "L2_Error_in_%: ",
            f"{l2normAbs * 100:.5f}%",
            "Rel_L2_Error_in_%: ",
            f"{l2normRel * 100:.5f}%",
        )
    else:
        print(
            "L2_Error: ",
            f"{l2normAbs:.5e}",
            " Rel_L2_Error: ",
            f"{l2normRel:.5e}",
        )


if __file__ == "__main__":
    cmdArgs = parseCommandLine()
    error = computeL2ErrorSquared(cmdArgs, **readL2ErrorData(cmdArgs))
    printL2Error(cmdArgs, **error)
