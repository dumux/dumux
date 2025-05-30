# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(test_ff_stokes_channel_3d EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_stokes_channel_3d PUBLIC "GRID_DIM=3")

dumux_add_test(NAME test_ff_stokes_channel_3d_cuboid
              TARGET test_ff_stokes_channel_3d
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_cuboid-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_cuboid-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d params.input
                             -Problem.Name test_ff_stokes_channel_3d_cuboid"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-12})

dumux_add_test(NAME test_ff_stokes_channel_3d_staircase
              TARGET test_ff_stokes_channel_3d
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK dune-subgrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_staircase-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_staircase-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d params.input
                             -Problem.Name test_ff_stokes_channel_3d_staircase -Problem.IsStaircaseGeometry true"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-12})

dumux_add_test(NAME test_ff_stokes_channel_3d_staircase_nocaching
               SOURCES main.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD HAVE_UMFPACK dune-subgrid_FOUND
               COMPILE_DEFINITIONS ENABLECACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_staircase-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_staircase_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d params.input
                             -Problem.Name test_ff_stokes_channel_3d_staircase_nocaching -Problem.IsStaircaseGeometry true"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-12})

dumux_add_test(NAME test_ff_stokes_channel_pseudo3d
              LABELS freeflow navierstokes
              SOURCES main.cc
              COMPILE_DEFINITIONS GRID_DIM=2
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_pseudo3d-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pseudo3d-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pseudo3d params_pseudo.input
                             -Problem.Name test_ff_stokes_channel_pseudo3d"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-12})


dune_symlink_to_source_files(FILES "params.input" "params_pseudo.input")
