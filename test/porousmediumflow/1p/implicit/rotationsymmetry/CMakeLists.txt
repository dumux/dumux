dune_symlink_to_source_files(FILES "params.input")

# using tpfa
dune_add_test(NAME test_1p_rotationsymmetry_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePRotSymTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_rotationsymmetry_tpfa-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_rotationsymmetry_tpfa-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_rotationsymmetry_tpfa params.input -Problem.Name test_1p_rotationsymmetry_tpfa")

# using box
dune_add_test(NAME test_1p_rotationsymmetry_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePRotSymBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_rotationsymmetry_box-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_rotationsymmetry_box-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_rotationsymmetry_box params.input -Problem.Name test_1p_rotationsymmetry_box")
