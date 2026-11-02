# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# non-isothermal tests
dumux_add_test(NAME test_solidenergy_tpfa
              LABELS porousmediumflow solidenergy
              SOURCES main.cc
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_solidenergy_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_solidenergy_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_solidenergy_tpfa params.input -Problem.Name test_solidenergy_tpfa")
