dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME test_projection_2d1d
              SOURCES test_projection_2d1d.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-functions_FOUND )"
              LABELS unit discretization
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzyData
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_tpfa_tpfa_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_tpfa_tpfa.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_tpfa_tpfa_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_tpfa_tpfa.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_box_box_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_box_box.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_box_box_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_box_box.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_fem_fem_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_fem_fem.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_fem_fem_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_fem_fem.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_tpfa_box_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_tpfa_box.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_tpfa_box_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_tpfa_box.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_box_tpfa_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_box_tpfa.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_box_tpfa_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_box_tpfa.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_fem_box_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_fem_box.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_fem_box_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_fem_box.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_fem_tpfa_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_fem_tpfa.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_fem_tpfa_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_fem_tpfa.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_tpfa_fem_2d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/2d_tpfa_fem.csv
                                ${CMAKE_SOURCE_DIR}/test/references/test_projection_2d1d_tpfa_fem_1d.csv
                                ${CMAKE_CURRENT_BINARY_DIR}/1d_tpfa_fem.csv
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_projection_2d1d params.input")
