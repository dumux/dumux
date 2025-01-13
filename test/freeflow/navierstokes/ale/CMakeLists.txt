# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

# Stokes version of the test (no inertia)
dumux_add_test(
        NAME test_ff_stokes_dfg_benchmark_movingmesh_pq1bubble_box
        SOURCES main.cc
        COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
        COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassBox
        LABELS freeflow navierstokes
        CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
        COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
        CMD_ARGS --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params.input
                     -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box -Problem.EnableInertiaTerms false
                     -LinearSolver.UseIterativeSolver true"
)
