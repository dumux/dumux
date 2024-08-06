# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
set(CMAKE_BUILD_TYPE Debug)
add_input_file_links()

dumux_add_test(NAME test_md_boundary_darcy2p2c_stokes1p2c_a02_a06_3d
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              SOURCES main.cc
              COMPILE_DEFINITIONS NONISOTHERMAL=0 DIM=3
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_stokes-00039.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_darcy-00039.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal params_3d.input
                             -Vtk.OutputName test_md_boundary_darcy2p2c_stokes1p2c_horizontal")

dumux_add_test(NAME test_md_boundary_darcy2p2cni_stokes1p2cni_a02_a06_3d
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              SOURCES main.cc
              COMPILE_DEFINITIONS NONISOTHERMAL=1 DIM=3
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_stokes-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_darcy-00040.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal params_nonisothermal_3d.input
                             -Vtk.OutputName test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal")


dumux_add_test(NAME test_md_boundary_darcy2p2c_stokes1p2c_a02_a06_2d
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              SOURCES main.cc
              COMPILE_DEFINITIONS NONISOTHERMAL=0 DIM=2
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_stokes-00039.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal_darcy-00039.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2c_stokes1p2c_horizontal params_2d.input
                             -Vtk.OutputName test_md_boundary_darcy2p2c_stokes1p2c_horizontal")

dumux_add_test(NAME test_md_boundary_darcy2p2cni_stokes1p2cni_a02_a06_2d
              LABELS multidomain multidomain_boundary stokesdarcy 2p2c navierstokes
              SOURCES main.cc
              COMPILE_DEFINITIONS NONISOTHERMAL=1 DIM=2
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_stokes-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_stokes-00040.vtu
                                     ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_darcy-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal_darcy-00040.vtu

                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal params_nonisothermal_2d.input
                             -Vtk.OutputName test_md_boundary_darcy2p2cni_stokes1p2cni_horizontal")
