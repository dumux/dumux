# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# comparison to analytical solution - only with tpfa
dumux_add_test(NAME test_richards_analytical_tpfa
              LABELS porousmediumflow richards
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsAnalyticalCC
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richards_analytical_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richards_analytical_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richards_analytical_tpfa params.input -Problem.Name test_richards_analytical_tpfa")
