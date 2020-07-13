dumux_add_test(NAME test_throattransmissibility
               SOURCES test_throattransmissibility.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_throattransmissibility)

dumux_add_test(NAME test_material_2p_porenetwork
               SOURCES test_material_2p_porenetwork.cc
               LABELS unit material
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzyData --delimiter " "
                        --files  ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_tetrahedron.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_tetrahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_octahedron.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_octahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_icosahedron.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_icosahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_dodecahedron.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_dodecahedron.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube_with_linear_regularization.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube_with_powerlaw_regularization.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube_with_powerlaw_regularization.dat
                                 ${CMAKE_SOURCE_DIR}/test/references/test_pnm_pcsw_cube_with_spline_regularization.dat
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_pcsw_cube_with_spline_regularization.dat
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_material_2p_porenetwork")
