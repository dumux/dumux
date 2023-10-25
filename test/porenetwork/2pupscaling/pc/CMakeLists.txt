# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(static)

add_input_file_links()

dumux_add_test(NAME test_pnm_2p_upscaling_pc_sw
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               COMPILE_DEFINITIONS CIRCLETHROAT=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00093.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p")
