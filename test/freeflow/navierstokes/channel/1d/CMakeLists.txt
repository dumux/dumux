# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_ff_navierstokes_1d
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_1d-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_1d-00001.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_1d params.input
                             -Problem.Name test_ff_navierstokes_1d")
dune_symlink_to_source_files(FILES params.input)
