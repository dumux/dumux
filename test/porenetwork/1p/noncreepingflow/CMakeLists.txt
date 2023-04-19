# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(nonisothermal)

add_input_file_links()

dumux_add_test(NAME test_pnm_1p_noncreeping_flow
              SOURCES main.cc
              LABELS porenetwork
              COMPILE_DEFINITIONS ISOTHERMAL=1
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p_noncreepingflow-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p_noncreepingflow-00000.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p_noncreeping_flow")
