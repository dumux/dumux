# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# if UG is found, we use UG for all tests. Otherwise, we use
# YaspGrid for the tests on unstructured grids and enforce
# the tests on unstructured grids to be skipped
if (dune-uggrid_FOUND)
    set(GRIDTYPE "Dune::UGGrid<2>")
else()
    set(GRIDTYPE "Dune::YaspGrid<2>")
endif()

dumux_add_test_target(NAME test_1p_convergence_analytic_tpfa
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=OnePConvergenceTpfa" "GRIDTYPE=${GRIDTYPE}")

dumux_add_test_target(NAME test_1p_convergence_analytic_mpfa
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=OnePConvergenceMpfa" "GRIDTYPE=${GRIDTYPE}")

dumux_add_test_target(NAME test_1p_convergence_analytic_box
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=OnePConvergenceBox" "GRIDTYPE=${GRIDTYPE}")

dumux_add_test_target(NAME test_1p_convergence_analytic_box_variables
                      SOURCES main.cc
                      COMPILE_DEFINITIONS "TYPETAG=OnePConvergenceBoxNewInterface" "GRIDTYPE=${GRIDTYPE}")

# using tpfa and structured grid with zero off-diagonal permeability entries
dumux_add_test(NAME test_1p_convergence_analytic_tpfa_structured
               TARGET test_1p_convergence_analytic_tpfa
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_tpfa params.input -Problem.C 0.0 -Problem.Name test_1p_convergence_analytic_tpfa_structured"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_analytic_tpfa_structured-errors.csv
                        --check "L2(p) mean 1.8")

# using mpfa and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_mpfa_structured
               TARGET test_1p_convergence_analytic_mpfa
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_mpfa params.input -Problem.Name test_1p_convergence_analytic_mpfa_structured"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_analytic_mpfa_structured-errors.csv
                        --check "L2(p) mean 1.8")

# using box and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_box_structured
               TARGET test_1p_convergence_analytic_box
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_box params.input -Problem.Name test_1p_convergence_analytic_box_structured"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_analytic_box_structured-errors.csv
                        --check "L2(p) mean 1.8")

# using box assembled with the variables per local dof and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_box_variables_structured
               TARGET test_1p_convergence_analytic_box_variables
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_box_variables params.input -Problem.Name test_1p_convergence_analytic_box_variables_structured"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_analytic_box_variables_structured-errors.csv
                        --check "L2(p) mean 1.8")

# tests on ustructured grids in case UG is found
if (dune-uggrid_FOUND)
    # using mpfa and unstructured grid
    dumux_add_test(NAME test_1p_convergence_analytic_mpfa_unstructured
                   TARGET test_1p_convergence_analytic_mpfa
                   LABELS porousmediumflow 1p
                   COMMAND ${DUMUX_CONVERGENCETEST}
                   CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_mpfa params.input -Grid.File ../../incompressible/grids/randomlydistorted.dgf -Problem.Name test_1p_convergence_analytic_mpfa_unstructured"
                            --levels 0 1 2 3
                            --level-args "-Grid.Refinement {}"
                            --errors test_1p_convergence_analytic_mpfa_unstructured-errors.csv
                            --check "L2(p) mean 1.8")

    # using box and unstructured grid
    dumux_add_test(NAME test_1p_convergence_analytic_box_unstructured
                   TARGET test_1p_convergence_analytic_box
                   LABELS porousmediumflow 1p
                   COMMAND ${DUMUX_CONVERGENCETEST}
                   CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_box params.input -Grid.File ../../incompressible/grids/randomlydistorted.dgf -Problem.Name test_1p_convergence_analytic_box_unstructured"
                            --levels 0 1 2 3
                            --level-args "-Grid.Refinement {}"
                            --errors test_1p_convergence_analytic_box_unstructured-errors.csv
                            --check "L2(p) mean 1.8")

# otherwise define test with the same name but enforce skip
else()
    dumux_add_test(NAME test_1p_convergence_analytic_mpfa_unstructured
                   SOURCES main.cc
                   CMAKE_GUARD dune-uggrid_FOUND
                   LABELS porousmediumflow 1p)

    dumux_add_test(NAME test_1p_convergence_analytic_box_unstructured
                   SOURCES main.cc
                   CMAKE_GUARD dune-uggrid_FOUND
                   LABELS porousmediumflow 1p)
endif()

# the schemes assembled with the new interfaces
dumux_add_test(NAME test_1p_convergence_analytic_pq3
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=OnePConvergencePQ3
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_pq3 params.input -Problem.Name test_1p_convergence_analytic_pq3"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_analytic_pq3-errors.csv
                        --check "L2(p) mean 3.6")

# PQ1Bubble is defined on simplices, so it needs a simplex grid
dumux_add_test(NAME test_1p_convergence_analytic_pq1bubble
               SOURCES main.cc
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS TYPETAG=OnePConvergencePQ1Bubble
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_pq1bubble params.input -Problem.Name test_1p_convergence_analytic_pq1bubble"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_analytic_pq1bubble-errors.csv
                        --check "L2(p) mean 1.8")

# the permeability on the faces between control volumes averaged instead of interpolated
dumux_add_test(NAME test_1p_convergence_analytic_pq1bubble_facetensoraverage
               TARGET test_1p_convergence_analytic_pq1bubble
               CMAKE_GUARD dune-alugrid_FOUND
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_analytic_pq1bubble params.input -Flux.ScvfPermeability FaceTensorAverage -Problem.Name test_1p_convergence_analytic_pq1bubble_facetensoraverage"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_analytic_pq1bubble_facetensoraverage-errors.csv
                        --check "L2(p) mean 1.8")
