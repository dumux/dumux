# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_pnm_2p_static
               SOURCES main.cc
               LABELS porenetwork
               COMMAND ${DUMUX_RUNTEST}
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzyData
                             --delimiter " "
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_static-reference.txt
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static_pc-s-curve.txt
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static")
