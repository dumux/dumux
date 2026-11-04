# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params_junction.input" "params_line.input")
dune_symlink_to_source_files(FILES "junction.dgf")

# with vertex control volumes on the median dual, the hybrid schemes of even order
# converge one order below the optimal rate in L2
dumux_add_test(NAME test_1p_convergence_line_pq2
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ2
               COMPILE_DEFINITIONS GRIDTYPE=Dune::OneDGrid
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_line_pq2 params_line.input -Problem.Name test_1p_convergence_line_pq2"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_line_pq2-errors.csv
                        --check "L2(p) finest 1.9 2.1"
                        --check "H1(p) finest 1.9 2.1")

dumux_add_test(NAME test_1p_convergence_line_pq3
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ3
               COMPILE_DEFINITIONS GRIDTYPE=Dune::OneDGrid
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_line_pq3 params_line.input -Problem.Name test_1p_convergence_line_pq3"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_line_pq3-errors.csv
                        --check "L2(p) finest 3.9 4.1"
                        --check "H1(p) finest 2.9 3.1")

dumux_add_test(NAME test_1p_convergence_junction_pq2
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ2
               COMPILE_DEFINITIONS GRIDTYPE=Dune::FoamGrid<1,3>
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_junction_pq2 params_junction.input -Problem.Name test_1p_convergence_junction_pq2"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_junction_pq2-errors.csv
                        --check "L2(p) finest 1.9 2.1"
                        --check "H1(p) finest 1.9 2.1")

dumux_add_test(NAME test_1p_convergence_junction_pq3
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ3
               COMPILE_DEFINITIONS GRIDTYPE=Dune::FoamGrid<1,3>
               LABELS porousmediumflow 1p
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_convergence_junction_pq3 params_junction.input -Problem.Name test_1p_convergence_junction_pq3"
                        --levels 0 1 2 3
                        --level-args "-Grid.Refinement {}"
                        --errors test_1p_convergence_junction_pq3-errors.csv
                        --check "L2(p) finest 3.9 4.1"
                        --check "H1(p) finest 2.9 3.1")
