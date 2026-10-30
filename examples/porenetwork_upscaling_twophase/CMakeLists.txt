# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME example_pnm2p_upscaling
               SOURCES main.cc
               LABELS porenetwork example
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzyData
                             --delimiter " "
                             --files ${CMAKE_SOURCE_DIR}/test/references/example_pnm2p_upscaling-reference.txt
                                     ${CMAKE_CURRENT_BINARY_DIR}/pnm_2p_upscaling_pc-s-curve.txt
                             --command "${CMAKE_CURRENT_BINARY_DIR}/example_pnm2p_upscaling params.input")
