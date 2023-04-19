# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(vertical)
dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_shallowwater_poiseuilleflow
               SOURCES main.cc
               LABELS shallowwater
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_poiseuilleflow-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/poiseuilleflow-00007.vtu
                              --zeroThreshold {"velocityY":1e-14}
                              --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_poiseuilleflow params.input")

dumux_add_test(NAME test_shallowwater_poiseuilleflow_parallel
               TARGET test_shallowwater_poiseuilleflow
               LABELS shallowwater
               CMAKE_GUARD MPI_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_poiseuilleflow-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/s0002-poiseuilleflow-parallel-00007.pvtu
                          --zeroThreshold {"velocityY":1e-14} --ignore "process rank"
                          --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_poiseuilleflow params.input -Problem.Name poiseuilleflow-parallel")

dumux_add_test(NAME test_shallowwater_poiseuilleflow_unstructured
               SOURCES main.cc
               LABELS shallowwater
               CMAKE_GUARD dune-uggrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_poiseuilleflow_unstructured-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/poiseuilleflow-unstructured-00007.vtu
                              --zeroThreshold {"velocityY":1e-14}
                              --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_poiseuilleflow_unstructured params.input -Problem.Name poiseuilleflow-unstructured -Grid.File grids/irregular_grid_10m.dgf")

dumux_add_test(NAME test_shallowwater_poiseuilleflow_unstructured_parallel
               TARGET test_shallowwater_poiseuilleflow_unstructured
               LABELS shallowwater
               CMAKE_GUARD "( MPI_FOUND AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_poiseuilleflow_unstructured-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/s0002-poiseuilleflow-unstructured-parallel-00007.pvtu
                          --zeroThreshold {"velocityY":1e-14} --ignore "process rank"
                          --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_poiseuilleflow_unstructured params.input -Problem.Name poiseuilleflow-unstructured-parallel -Grid.File grids/irregular_grid_10m.dgf")

dune_symlink_to_source_files(FILES "grids")
