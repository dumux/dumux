# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES grids "params.input")

dumux_add_test(NAME test_md_dualnetwork_heat_1p_solid
              LABELS multidomain dualnetwork 1p
              SOURCES main.cc
              CMAKE_GUARD "( dune-subgrid_FOUND AND dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_dualnetwork_heat_solid-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_dualnetwork_heat_solid-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_dualnetwork_heat_void-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_dualnetwork_heat_void-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_dualnetwork_heat_1p_solid params.input
                                    -Vtk.OutputName test_md_dualnetwork_heat")
