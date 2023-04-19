# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c
              LABELS multidomain multidomain_embedded 1pnc richardsnc
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-uggrid_FOUND )"
              TIMEOUT 1500
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p2c_richards2c_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p2c_richards2c_1d-00003.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p2c_richards2c_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p2c_richards2c_3d-00003.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p2c_richards2c params.input
                        -Vtk.OutputName test_md_embedded_1d3d_1p2c_richards2c")
dune_symlink_to_source_files(FILES "params.input")
