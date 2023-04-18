# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_ff_navierstokes_angeli
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_angeli-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_angeli-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_angeli params.input
                             -Problem.Name test_ff_navierstokes_angeli")

dumux_add_test(NAME test_ff_navierstokes_angeli_averaged
              LABELS freeflow navierstokes
              TARGET test_ff_navierstokes_angeli
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_angeli_averaged-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_angeli_averaged-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_angeli params.input
                             -TimeLoop.TEnd 1e-6 -TimeLoop.DtInitial 1e-6
                             -Problem.Name test_ff_navierstokes_angeli_averaged
                             -Problem.InterpolateExactVelocity true")

dune_symlink_to_source_files(FILES "params.input")
