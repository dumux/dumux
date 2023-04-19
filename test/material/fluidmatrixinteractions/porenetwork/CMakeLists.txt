# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_throattransmissibility
               LABELS unit material porenetwork
               SOURCES test_throattransmissibility.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_throattransmissibility)

dumux_add_test(NAME test_material_2p_porenetwork
               SOURCES test_material_2p_porenetwork.cc
               LABELS unit material porenetwork
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData --delimiter " "
                        --files  ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_tetrahedron-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_tetrahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_octahedron-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_octahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_icosahedron-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_icosahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_dodecahedron-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_dodecahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube_with_linear_regularization.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube_with_powerlaw_regularization-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube_with_powerlaw_regularization.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube_with_spline_regularization-reference.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube_with_spline_regularization.dat
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_material_2p_porenetwork")
