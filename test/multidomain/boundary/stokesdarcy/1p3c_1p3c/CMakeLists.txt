# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_md_boundary_darcy1p3c_stokes1p3c_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 1pnc navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p3c_stokes1p3c_horizontal_stokes-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p3c_stokes1p3c_horizontal_stokes-00015.vtu
                                    ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p3c_stokes1p3c_horizontal_darcy-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p3c_stokes1p3c_horizontal_darcy-00015.vtu

                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p3c_stokes1p3c_horizontal params.input")
