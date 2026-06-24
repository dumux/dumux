# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
set(CMAKE_BUILD_TYPE Debug)
add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME example_pnm2p_upscaling_drainage
               SOURCES main.cc
               COMPILE_DEFINITIONS DRAINAGE=1
               LABELS porenetwork
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzyData
                             --delimiter " "
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_static-reference.txt
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static_pc-s-curve.txt
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static")


dumux_add_test(NAME example_pnm2p_upscaling_imbibition
               SOURCES main.cc
               COMPILE_DEFINITIONS DRAINAGE=0
               LABELS porenetwork
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzyData
                       --delimiter " "
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_static-reference.txt
                               ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static_pc-s-curve.txt
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static")
