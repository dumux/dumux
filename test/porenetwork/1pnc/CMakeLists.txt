# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

add_executable(test_pnm_1p2c EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_pnm_1p2c PRIVATE ISOTHERMAL=1)

dumux_add_test(NAME test_pnm_1p2c_advection
               TARGET test_pnm_1p2c
               LABELS porenetwork
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p2c-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c-00055.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c")

dumux_add_test(NAME test_pnm_1p2c_diffusion
               TARGET test_pnm_1p2c
               LABELS porenetwork
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p2c_diffusion-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c_diffusion-00054.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2c params_diffusion.input")

dumux_add_test(NAME test_pnm_1p2cni
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_1p2cni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2cni-00056.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_1p2cni -Problem.Name test_pnm_1p2cni")
