# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_2p_buckleyleverett_tpfa
               LABELS porousmediumflow 2p
               SOURCES main.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_2p_buckleyleverett_tpfa
               CMD_ARGS params.input -Problem.Name test_2p_buckleyleverett_tpfa)
