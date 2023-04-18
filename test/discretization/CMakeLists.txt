# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(cellcentered)
add_subdirectory(facecentered)
add_subdirectory(staggered)
add_subdirectory(box)
add_subdirectory(pq1bubble)
add_subdirectory(projection)
add_subdirectory(rotationsymmetry)
add_input_file_links()

add_executable(test_walldistance EXCLUDE_FROM_ALL test_walldistance.cc)

dumux_add_test(NAME test_walldistance_2dcube
               TARGET test_walldistance
               LABELS unit
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "wallElementIndex"
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_walldistance params_walldistance.input -TestCase 2DCube"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2dcube_box_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_2dcube_box_element.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2dcube_box_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_2dcube_box_vertex.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2dcube_cctpfa_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_2dcube_cctpfa_element.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2dcube_cctpfa_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_2dcube_cctpfa_vertex.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_walldistance_3dcube
               TARGET test_walldistance
               LABELS unit
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "wallElementIndex"
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_walldistance params_walldistance.input -TestCase 3DCube"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_box_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dcube_box_element.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_box_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dcube_box_vertex.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_cctpfa_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dcube_cctpfa_element.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_cctpfa_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dcube_cctpfa_vertex.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_walldistance_3dcube_parallel
               TARGET test_walldistance
               LABELS unit
               CMAKE_GUARD MPI_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "wallElementIndex"
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_walldistance params_walldistance.input -TestCase 3DCube"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_box_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_aabb_tree_3dcube_box_element.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_box_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_aabb_tree_3dcube_box_vertex.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_cctpfa_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_aabb_tree_3dcube_cctpfa_element.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dcube_cctpfa_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_aabb_tree_3dcube_cctpfa_vertex.pvtu
                        --relative 1e-6)

dune_symlink_to_source_files(FILES heart.msh) # grid file
dumux_add_test(NAME test_walldistance_3dheart
               TARGET test_walldistance
               LABELS unit
               CMAKE_GUARD dune-alugrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "wallElementIndex"
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_walldistance params_walldistance.input -TestCase 3DMesh"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dmesh_box_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dmesh_box_element.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dmesh_box_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dmesh_box_vertex.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dmesh_cctpfa_element-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dmesh_cctpfa_element.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3dmesh_cctpfa_vertex-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3dmesh_cctpfa_vertex.vtu
                        --relative 1e-6)
