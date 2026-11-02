# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_1p_pointsources_timedependent_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMMAND ${DUMUX_RUNTEST}
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityTimeDependentCCTpfa
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timedependent_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timedependent_tpfa-00003.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timedependent_tpfa params.input -TimeLoop.TEnd 4
                       -Problem.Name test_1p_pointsources_timedependent_tpfa")
