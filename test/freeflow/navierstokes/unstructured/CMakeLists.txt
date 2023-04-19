# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input cylinder_channel.msh cylinder_channel_quad.msh)

# Navier-Stokes version of the test (Re=20)
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumDiamond
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassTpfa
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_stationary_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa -Problem.EnableInertiaTerms true")

# Stokes version of the test (no inertia)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_diamond_tpfa
               TARGET test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_diamond_tpfa-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_diamond_tpfa -Problem.EnableInertiaTerms false
                             -LinearSolver.UseIterativeSolver true")

# Navier-Stokes version of the test (Re=20)
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassBox
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box -Problem.EnableInertiaTerms true")

# Stokes version of the test (no inertia)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box
               TARGET test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box -Problem.EnableInertiaTerms false
                             -LinearSolver.UseIterativeSolver true")

# Navier-Stokes version of the test (Re=20)
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassDiamond
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond -Problem.EnableInertiaTerms true
                             -Grid.File cylinder_channel_quad.msh")

# Stokes version of the test (no inertia)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond
               TARGET test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond -Problem.EnableInertiaTerms false
                             -Grid.File cylinder_channel_quad.msh")
