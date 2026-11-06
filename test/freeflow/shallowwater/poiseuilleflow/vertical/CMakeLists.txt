# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_shallowwater_poiseuilleflow_vertical
               SOURCES main.cc
               LABELS shallowwater
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_shallowwater_poiseuilleflow_vertical-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/poiseuilleflow-00001.vtu
                              --zeroThreshold {"velocityY":1e-14}
                              --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_poiseuilleflow_vertical params.input")

# convergence towards the exact mean velocity of the film flow
dumux_add_test(NAME test_shallowwater_poiseuilleflow_vertical_convergence
               TARGET test_shallowwater_poiseuilleflow_vertical
               LABELS shallowwater
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_poiseuilleflow_vertical params.input -Problem.Name test_shallowwater_poiseuilleflow_vertical_convergence"
                        --levels 20 40 80
                        --level-args "-Grid.Cells '{} 2'"
                        --errors test_shallowwater_poiseuilleflow_vertical_convergence-errors.csv
                        --check "L2Rel(u) mean 0.9")
