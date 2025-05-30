# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_ff_navierstokes_kovasznay
              LABELS freeflow navierstokes
              TIMEOUT 3600
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS UPWINDSCHEMEORDER=1
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_kovasznay-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay params.input
                             -Problem.Name test_ff_navierstokes_kovasznay")

add_executable(test_ff_navierstokes_kovasznay_higherorder EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_navierstokes_kovasznay_higherorder PUBLIC "UPWINDSCHEMEORDER=2")

dumux_add_test(NAME test_ff_navierstokes_kovasznay_higherorder_cuboid
              TARGET test_ff_navierstokes_kovasznay_higherorder
              LABELS freeflow navierstokes
              TIMEOUT 3600
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_kovasznay_higherorder_cuboid-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay_higherorder_cuboid-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay_higherorder params.input
                             -Problem.Name test_ff_navierstokes_kovasznay_higherorder_cuboid
                             -Flux.TvdApproach Hou
                             -Flux.DifferencingScheme Vanleer")

dumux_add_test(NAME test_ff_navierstokes_kovasznay_higherorder_staircase
              TARGET test_ff_navierstokes_kovasznay_higherorder
              LABELS freeflow navierstokes
              TIMEOUT 3600
              CMAKE_GUARD "( HAVE_UMFPACK AND dune-subgrid_FOUND )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_kovasznay_higherorder_staircase-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay_higherorder_staircase-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay_higherorder params.input
                             -Problem.Name test_ff_navierstokes_kovasznay_higherorder_staircase
                             -Problem.IsStaircaseGeometry true
                             -Flux.TvdApproach Hou
                             -Flux.DifferencingScheme Vanleer")

dune_symlink_to_source_files(FILES "params.input")
