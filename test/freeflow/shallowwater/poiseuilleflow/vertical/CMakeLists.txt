# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_shallowwater_poiseuilleflow_vertical
               SOURCES main.cc
               LABELS shallowwater
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_shallowwater_poiseuilleflow_vertical-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/poiseuilleflow-00001.vtu
                              --zeroThreshold {"velocityY":1e-14}
                              --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_poiseuilleflow_vertical params.input")
