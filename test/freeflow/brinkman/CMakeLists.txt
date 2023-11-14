# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_ff_brinkman_obstacle
              LABELS freeflow navierstokes brinkman
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_brinkman_obstacle-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_brinkman_obstacle-00000.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_brinkman_obstacle params.input
                             -Problem.Name test_ff_brinkman_obstacle")
