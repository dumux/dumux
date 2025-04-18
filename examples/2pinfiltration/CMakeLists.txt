# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "initialsolutioncc.txt")

# using tpfa and point source
dumux_add_test(NAME example_2p_pointsource_adaptive
               LABELS porousmediumflow 2p example
               SOURCES main.cc
               CMAKE_GUARD dune-alugrid_FOUND
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy --zeroThreshold {"pc":0.1}
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_2padaptive-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/example_2p_pointsource_adaptive-00006.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_2p_pointsource_adaptive params.input -TimeLoop.TEnd 100 -TimeLoop.DtInitial 10 -Problem.Name  example_2p_pointsource_adaptive")
