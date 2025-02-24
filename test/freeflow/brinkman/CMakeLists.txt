# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_ff_brinkman_obstacle_staggered_tpfa
               LABELS freeflow navierstokes brinkman
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_brinkman_obstacle_staggered_tpfa-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_brinkman_obstacle_staggered_tpfa-00000.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_brinkman_obstacle_staggered_tpfa params.input
                             -Problem.Name test_ff_brinkman_obstacle_staggered_tpfa")

dumux_add_test(NAME test_ff_brinkman_obstacle_pq1bubble_box_simplex
               LABELS freeflow navierstokes brinkman
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ1BubbleModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_brinkman_obstacle_pq1bubble_box_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_brinkman_obstacle_pq1bubble_box_simplex-00000.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_brinkman_obstacle_pq1bubble_box_simplex params.input
                             -Problem.Name test_ff_brinkman_obstacle_pq1bubble_box_simplex")
