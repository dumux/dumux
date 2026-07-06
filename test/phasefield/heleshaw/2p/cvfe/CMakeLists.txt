# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_heleshaw_2p_cvfe_saffman_taylor
               LABELS phasefield heleshaw
               SOURCES main.cc
               CMAKE_GUARD DUMUX_HAVE_GRIDFORMAT
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_heleshaw_2p_cvfe_saffman_taylor-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_heleshaw_2p_cvfe_saffman_taylor-00021.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_heleshaw_2p_cvfe_saffman_taylor params.input")
