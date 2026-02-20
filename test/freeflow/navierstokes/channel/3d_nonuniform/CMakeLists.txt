# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "channel.msh" "channel_quad.msh")
dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform_diamond
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=ThreeDChannelTestMomentumDiamond
               COMPILE_DEFINITIONS TYPETAG_MASS=ThreeDChannelTestMassTpfa
               LABELS freeflow navierstokes
               TIMEOUT 3600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform_diamond-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_diamond-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_diamond params.input
                             -Problem.Name test_ff_stokes_channel_3d_nonuniform_diamond")

dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform_diamond_simplex
               TARGET test_ff_stokes_channel_3d_nonuniform_diamond
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform_diamond_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_diamond_simplex-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_diamond params.input
                             -Problem.Name test_ff_stokes_channel_3d_nonuniform_diamond_simplex
                             -Grid.File channel.msh")

dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform_pq1bubble_simplex
                SOURCES main.cc
                COMPILE_DEFINITIONS TYPETAG_MOMENTUM=ThreeDChannelTestMomentumPQ1Bubble
                COMPILE_DEFINITIONS TYPETAG_MASS=ThreeDChannelTestMassBox
                COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
                COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
                COMPILE_DEFINITIONS USE_STOKES_SOLVER=1
                LABELS freeflow navierstokes
                CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
                COMMAND ${DUMUX_RUNTEST}
                CMD_ARGS      --script fuzzy
                                --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform_pq1bubble_simplex-reference.vtu
                                        ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_pq1bubble_simplex-00001.vtu
                                --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_pq1bubble_simplex params.input
                                -Problem.Name test_ff_stokes_channel_3d_nonuniform_pq1bubble_simplex
                                -Grid.File channel.msh -Assembly.Multithreading false")

dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform_pq2
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=ThreeDChannelTestMomentumPQ2
               COMPILE_DEFINITIONS TYPETAG_MASS=ThreeDChannelTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               LABELS freeflow navierstokes
               TIMEOUT 600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform_pq2-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_pq2-00001.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_pq2 params.input
                               -Problem.Name test_ff_stokes_channel_3d_nonuniform_pq2")

dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform_pq2_simplex
               TARGET test_ff_stokes_channel_3d_nonuniform_pq2
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform_pq2_simplex-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_pq2_simplex-00001.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_pq2 params.input
                               -Problem.Name test_ff_stokes_channel_3d_nonuniform_pq2_simplex
                               -Grid.File channel.msh")
