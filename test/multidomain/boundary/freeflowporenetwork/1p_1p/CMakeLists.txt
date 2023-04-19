# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_md_boundary_ff1p_pnm1p
              LABELS multidomain multidomain_boundary freeflowpnm 1p navierstokes porenetwork
              SOURCES main.cc
              CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p_pnm1p_ff-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pnm1p_ff-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p_pnm1p_pnm.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pnm1p_pnm-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pnm1p params.input
                                   -Vtk.OutputName test_md_boundary_ff1p_pnm1p")

dumux_add_test(NAME test_md_boundary_ff1p_pnm1p_single_throat
              LABELS multidomain multidomain_boundary freeflowpnm 1p navierstokes porenetwork
              TARGET test_md_boundary_ff1p_pnm1p
              CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p_pnm1p_single_throat_ff-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pnm1p_single_throat_ff-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1p_pnm1p_single_throat_pnm-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pnm1p_single_throat_pnm-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pnm1p params_single_throat.input
                                   -Vtk.OutputName test_md_boundary_ff1p_pnm1p_single_throat")
