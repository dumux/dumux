#!/usr/bin/env python3

import os
import sys

from argparse import ArgumentParser

try:
    commonPath = os.path.dirname(os.path.abspath(__file__))
    # It is assumed that functions are provided by file in parent folder
    commonPath = commonPath.rsplit("/build-cmake")[0] + commonPath.rsplit("/build-cmake")[1] + "/.."
    sys.path.append(commonPath)
    from convergencetests import runTestsAndPlotResults
except Exception:
    sys.exit("Could not import runTestsAndPlotResults")


#sys.path.append('/home/martins/DUMUX/dumux/test/freeflow/navierstokes')


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Run script for the convergence test."
    )
    parser.add_argument("-n", "--num-refinements", required=False, default=5)
    #parser.add_argument("-e", "--executable", type=str, required=True)
    #parser.add_argument("-o", "--outputname", type=str, required=True)
    args = vars(parser.parse_args())

    numRefinements = int(args["num_refinements"])

    testPrefix = "donea-momentum"
    executables = ["test_ff_stokes_donea_momentum_diamond_simplex", "test_ff_stokes_donea_momentum_diamond_quad"]
    grids =    [
                [ ["simplex-structured", []], ["simplex-unstructured", ["-Grid.File ./grids/unstructured_simplex.msh"]] ],
                [ ["quad-structured", []],    ["quad-unstructured",    ["-Grid.File ./grids/unstructured_quad.msh"   ]] ]
               ]

    methods =    [
                  ["weak-sym", []],
                  ["unsym", ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]]
                 ]

    runTestsAndPlotResults(testPrefix, executables, grids, methods, numRefinements)
