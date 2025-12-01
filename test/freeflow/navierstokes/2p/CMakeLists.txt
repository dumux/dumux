# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_ff_stokes_2p_pq1bubble_simplex
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=ThreeDChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=ThreeDChannelTestMassBox
               LABELS freeflow navierstokes
               TIMEOUT 3600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_2p_pq1bubble_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_2p_pq1bubble_simplex-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_2p_pq1bubble_simplex params.input
                             -Problem.Name test_ff_stokes_2p_pq1bubble_simplex")
