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

    testNames = ["sincos-momentum-simplex-structured",
                 "sincos-momentum-simplex-unstructured",
                 "sincos-momentum-quad-structured",
                 "sincos-momentum-quad-unstructured"]

    testRuns = [
                [ ["test_ff_navierstokes_sincos_momentum_diamond_simplex", []] ],
                [ ["test_ff_navierstokes_sincos_momentum_diamond_simplex", ["-Grid.File ../donea/grids/unstructured_simplex.msh"]] ],
                [ ["test_ff_navierstokes_sincos_momentum_diamond_quad", []], ["test_ff_navierstokes_sincos_momentum_staggered_quad", []] ],
                [ ["test_ff_navierstokes_sincos_momentum_diamond_quad", ["-Grid.File ../donea/grids/unstructured_quad.msh"   ]] ]
               ]

    subRuns = [
                [ [ ["weak-sym", []], ["unsym", ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ] ],
                [ [ ["weak-sym", []], ["unsym", ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ] ],
                [ [ ["weak-sym", []], ["unsym", ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ], [ ["staggered", []] ] ],
                [ [ ["weak-sym", []], ["unsym", ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ] ]
              ]

    runTestsAndPlotResults(testNames, "params_momentum.input", testRuns, subRuns, numRefinements)
