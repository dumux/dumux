# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_phasefield_conservativelevelset_pq2
               LABELS phasefield conservativelevelset
               SOURCES main.cc
               CMAKE_GUARD DUMUX_HAVE_GRIDFORMAT
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_phasefield_conservativelevelset_pq2-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_phasefield_conservativelevelset_pq2-00023.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_phasefield_conservativelevelset_pq2 params.input
                        -Grid.Cells \" 10 10\" -Output.Every 1000")
