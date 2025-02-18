# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

########## ISOTHERMAL ##########
dumux_add_test(NAME test_md_boundary_ff1p2c_pnm1p2c
                LABELS multidomain multidomain_boundary freeflowpnm 1pnc navierstokes porenetwork
                SOURCES main.cc
                COMPILE_DEFINITIONS ISOTHERMAL=1
                CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMD_ARGS  --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p2c_pnm1p2c_ff-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2c_pnm1p2c_ff-00019.vtu
                                  ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p2c_pnm1p2c_pnm-reference.vtp
                                  ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2c_pnm1p2c_pnm-00019.vtp
                          --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2c_pnm1p2c params.input
                                     -Vtk.OutputName test_md_boundary_ff1p2c_pnm1p2c")

dumux_add_test(NAME test_md_boundary_ff1p2c_pnm1p2c_diffusion
                LABELS multidomain multidomain_boundary freeflowpnm 1pnc navierstokes porenetwork
                SOURCES main.cc
                COMPILE_DEFINITIONS ISOTHERMAL=1
                CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMD_ARGS  --script fuzzy
                          --ignore "velocity_gas (m/s)"
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p2c_pnm1p2c_diffusion_ff-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2c_pnm1p2c_diffusion_ff-00018.vtu
                                  ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p2c_pnm1p2c_diffusion_pnm-reference.vtp
                                  ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2c_pnm1p2c_diffusion_pnm-00018.vtp
                          --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2c_pnm1p2c_diffusion params.input
                                    -Problem.OnlyDiffusion true -PNM.Problem.InitialPressure 1e5 -Vtk.OutputName test_md_boundary_ff1p2c_pnm1p2c_diffusion")

########## NON-ISOTHERMAL #########
dumux_add_test(NAME test_md_boundary_ff1p2cni_pnm1p2cni
                LABELS multidomain multidomain_boundary freeflowpnm 1pnc navierstokes porenetwork
                SOURCES main.cc
                COMPILE_DEFINITIONS ISOTHERMAL=0
                CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
                COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p2cni_pnm1p2cni_ff-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2cni_pnm1p2cni_ff-00019.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p2cni_pnm1p2cni_pnm-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2cni_pnm1p2cni_pnm-00019.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p2cni_pnm1p2cni paramsni.input
                                    -Vtk.OutputName test_md_boundary_ff1p2cni_pnm1p2cni")
