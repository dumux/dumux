# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "params_projection.input" "tinyroot.dgf" "tinyroot_3d.msh")

dumux_add_test(NAME test_md_embedded_1d3d_1p_richards_tpfatpfa
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfatpfa_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfatpfa_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfatpfa_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfatpfa params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p_richards_tpfatpfa")

dumux_add_test(NAME test_md_embedded_1d3d_1p_richards_tpfabox
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-foamgrid_FOUND
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfabox_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_tpfabox_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfabox_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_tpfabox params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p_richards_tpfabox")

# Test for the projection coupling manager (a resolved root soil interface 1D-3D method)
dumux_add_test(NAME test_md_embedded_1d3d_1p_richards_proj
              LABELS multidomain multidomain_embedded 1p richards
              SOURCES main.cc
              COMPILE_DEFINITIONS SOILTYPETAG=SoilBox COUPLINGMODE=Projection SOILGRID=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TIMEOUT 1500
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_proj_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_proj_1d-00004.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p_richards_proj_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_proj_3d-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p_richards_proj params_projection.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p_richards_proj")
