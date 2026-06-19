# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_heleshaw_2p_saffman_taylor
               SOURCES main.cc
               LABELS freeflow heleshaw
               TIMEOUT 600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_heleshaw_2p_saffman_taylor
                                   params.input -TimeLoop.TEnd 0.05
                                   -Problem.Name test_heleshaw_2p_saffman_taylor")
