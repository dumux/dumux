# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(convergence)
add_subdirectory(transient)

add_input_file_links()

add_executable(test_md_boundary_ff1p_pm1p EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_horizontal
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              TARGET test_md_boundary_ff1p_pm1p
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_horizontal_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_horizontal_freeflow-00000.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_horizontal_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_horizontal_darcy-00000.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p params_horizontalflow.input
                                   -Vtk.OutputName test_md_boundary_ff1p_pm1p_horizontal")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_vertical
              LABELS multidomain multidomain_boundary stokesdarcy 1p navierstokes
              TARGET test_md_boundary_ff1p_pm1p
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_vertical_stokes-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_vertical_freeflow-00000.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_stokes1p_vertical_darcy-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_vertical_darcy-00000.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p params_verticalflow.input
                                   -Vtk.OutputName test_md_boundary_ff1p_pm1p_vertical"
                        --zeroThreshold {"velocity_liq \(m/s\)_0":6e-17})
