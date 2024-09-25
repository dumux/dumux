# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(static)

add_input_file_links()

dumux_add_test(NAME test_pnm_2p
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p-00171.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p")

dumux_add_test(NAME test_pnm_2p_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 USETHETAREGULARIZATION=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_reg-reference.vtp
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_reg-00259.vtp
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_reg -Problem.Name test_pnm_2p_reg")

dumux_add_test(NAME test_pnm_2pni
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni-00133.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pni params_ni.input -Problem.Name test_pnm_2pni")
