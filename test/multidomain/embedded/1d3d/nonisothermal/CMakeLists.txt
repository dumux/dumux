# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# line coupling mode
dumux_add_test(NAME test_wellbore_heat_transport_line
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=SoilCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_bm_singlepipe_surface_tpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_bm_singlepipe_surface_tpfa_1d-00002.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_bm_singlepipe_surface_tpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_bm_singlepipe_surface_tpfa_3d-00002.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_benchmark_surface_tpfa singlepipeparams.input -TimeLoop.TEnd 31536000 -Vtk.OutputName test_bm_singlepipe_surface_tpfa")


# surface coupling mode
dumux_add_test(NAME test_wellbore_heat_transport_surface
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=SoilCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/dumux_runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_bm_singlepipe_surface_tpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_bm_singlepipe_surface_tpfa_1d-00002.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_bm_singlepipe_surface_tpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_bm_singlepipe_surface_tpfa_3d-00002.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_benchmark_surface_tpfa singlepipeparams.input -TimeLoop.TEnd 31536000 -Vtk.OutputName test_bm_singlepipe_surface_tpfa")


add_input_file_links()
dune_symlink_to_source_files(FILES "grids")
