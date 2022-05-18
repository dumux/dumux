#!/usr/bin/env python3

import os
import sys

from argparse import ArgumentParser

try:
    commonPath = os.path.dirname(os.path.abspath(__file__))
    # It is assumed that functions are provided by file in parent folder
    commonPath = commonPath.rsplit("/build-cmake")[0] + commonPath.rsplit("/build-cmake")[1] + "/.."
    sys.path.append(commonPath)
    from convergencetests import runTests
    from convergencetests import plotResults
except Exception:
    sys.exit("Could not import functions for running tests")


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

    def refinement(ref):
        return " -Grid.Refinement " + str(ref)

    def delaunay(ref):
        return "-Grid.File ../donea/grids/delaunay_" + str(ref) + ".msh"

    def unstructuredSimplex(ref):
        return "-Grid.File ../donea/grids/unstructured_simplex.msh " + refinement(ref)

    def unstructuredQuad(ref):
        return "-Grid.File ../donea/grids/unstructured_quad.msh " + refinement(ref)

    ########### only momentum balance ###########
    #############################################

    # testNames = ["donea-momentum-simplex-delaunay",
    #              "donea-momentum-simplex-unstructured",
    #              "donea-momentum-quad-structured",
    #              "donea-momentum-quad-unstructured"]

    # testRuns = [
    #             [ ["test_ff_stokes_donea_momentum_diamond_simplex", delaunay] ],
    #             [ ["test_ff_stokes_donea_momentum_diamond_simplex", unstructuredSimplex] ],
    #             [ ["test_ff_stokes_donea_momentum_diamond_quad", refinement], ["test_ff_stokes_donea_momentum_staggered_quad", refinement] ],
    #             [ ["test_ff_stokes_donea_momentum_diamond_quad", unstructuredQuad] ]
    #            ]

    # subRuns = [
    #             [ [ [["donea_diamond_momentum_simplex_struct",        "weak-sym"], []], [["donea_diamond_momentum_simplex_struct",   "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ] ],
    #             [ [ [["donea_diamond_momentum_simplex_unstruct",      "weak-sym"], []], [["donea_diamond_momentum_simplex_unstruct", "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ] ],
    #             [ [ [["donea_diamond_momentum_quad_struct",           "weak-sym"], []], [["donea_diamond_momentum_quad_struct",      "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ], [ [["donea_staggered_momentum_quad_struct",           "staggered"], []] ] ],
    #             [ [ [["donea_diamond_momentum_quad_unstruct",         "weak-sym"], []], [["donea_diamond_momentum_quad_unstruct",    "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient", "true"]] ] ]
    #           ]

    # runTests(testNames, "params.input", testRuns, subRuns, numRefinements)
    # plotResults(testNames, testRuns, subRuns, [ [[2,3],"v"] ])

    ########### full navier-stokes ###########
    ##########################################

    method = "tpfa"
    testNames = ["donea-simplex-delaunay",
                 "donea-simplex-unstructured",
                 "donea-quad-structured",
                 "donea-quad-unstructured"]

    testRuns = [
                [ ["test_ff_stokes_donea_diamond_"+method+"_simplex", delaunay] ],
                [ ["test_ff_stokes_donea_diamond_"+method+"_simplex", unstructuredSimplex] ],
                [ ["test_ff_stokes_donea_diamond_"+method+"_quad", refinement], ["test_ff_stokes_donea", refinement] ],
                [ ["test_ff_stokes_donea_diamond_"+method+"_quad", unstructuredQuad] ]
               ]

    subRuns = [
                [ [ [["donea_diamond_"+method+"_simplex_struct",   "weak-sym"], []], [["donea_diamond_"+method+"_simplex_struct",   "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient true"]] ] ],
                [ [ [["donea_diamond_"+method+"_simplex_unstruct", "weak-sym"], []], [["donea_diamond_"+method+"_simplex_unstruct", "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient true"]] ] ],
                [ [ [["donea_diamond_"+method+"_quad_struct",      "weak-sym"], []], [["donea_diamond_"+method+"_quad_struct",      "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient true"]] ], [ [["donea_staggered_tpfa_quad_struct",   "staggered"], []] ] ],
                [ [ [["donea_diamond_"+method+"_quad_unstruct",    "weak-sym"], []], [["donea_diamond_"+method+"_quad_unstruct",    "unsym"], ["-FreeFlow.EnableUnsymmetrizedVelocityGradient true"]] ] ]
              ]

    runTests(testNames, "params.input", testRuns, subRuns, numRefinements)
    plotResults(testNames, testRuns, subRuns, [[[2,5],"p",1], [[4,6],"v",2] ])