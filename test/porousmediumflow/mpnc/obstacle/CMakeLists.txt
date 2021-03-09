add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_mpnc_obstacle_box
              LABELS porousmediumflow mpnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ObstacleBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_obstacle_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_obstacle_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_obstacle_box params.input -Problem.Name test_mpnc_obstacle_box")

dumux_add_test(NAME test_mpnc_obstacle_tpfa
              LABELS porousmediumflow mpnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ObstacleCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_obstacle_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_obstacle_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_obstacle_tpfa params.input -Problem.Name test_mpnc_obstacle_tpfa")
