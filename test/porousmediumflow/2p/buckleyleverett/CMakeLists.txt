# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_2p_buckleyleverett_tpfa
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_buckleyleverett_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_buckleyleverett_tpfa-00035.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_buckleyleverett_tpfa params.input -Problem.Name test_2p_buckleyleverett_tpfa")
