add_subdirectory(cellcentered)
add_subdirectory(facecentered)
add_subdirectory(staggered)
add_subdirectory(box)
add_subdirectory(projection)
add_subdirectory(rotationsymmetry)
add_input_file_links()

dumux_add_test(NAME test_disc_fvgridvariables
               LABELS unit
               SOURCES test_fvgridvariables.cc
               COMMAND ./test_disc_fvgridvariables
               CMD_ARGS -Problem.Name gridvarstest)

add_executable(test_walldistance EXCLUDE_FROM_ALL test_walldistance.cc)

dumux_add_test(NAME test_walldistance
               TARGET test_walldistance
               LABELS unit
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_walldistance params_walldistance.input"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_bounding_spheres_2D_cctpfa.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_brute_force_2D_cctpfa.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_2D_cctpfa.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_bounding_spheres_2D_box.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_brute_force_2D_box.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_2D_box.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_bounding_spheres_3D_cctpfa.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_brute_force_3D_cctpfa.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3D_cctpfa.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_bounding_spheres_3D_box.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_brute_force_3D_box.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_aabb_tree_3D_box.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_walldistance_parallel
               TARGET test_walldistance
               LABELS unit
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --zeroThreshold {"wallElementIndex":1000}
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_walldistance params_walldistance.input"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_bounding_spheres_2D_cctpfa.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_2D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_bounding_spheres_2D_box.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_bounding_spheres_3D_cctpfa.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_3D_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-result_bounding_spheres_3D_box.pvtu
                        --relative 1e-6)

dumux_add_test(NAME test_walldistance_cake
               TARGET test_walldistance
               LABELS unit
               CMAKE_GUARD "( dune-alugrid_FOUND OR dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_walldistance params_walldistance.input -Test.RunCakeGrid true"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_cake-cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_bounding_spheres_3D_cctpfa_cake.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_cake-cctpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_brute_force_3D_cctpfa_cake.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_cake-box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_bounding_spheres_3D_box_cake.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_walldistance_cake-box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/result_brute_force_3D_box_cake.vtu
                        --relative 1e-6)
