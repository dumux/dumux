add_input_file_links()

add_executable(test_1pni_simple_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1pni_simple_tpfa PUBLIC TYPETAG=OnePNISimpleTpfa)

add_executable(test_1pni_simple_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1pni_simple_box PUBLIC TYPETAG=OnePNISimpleBox)

dumux_add_test(NAME test_1pni_simple_linear_tpfa
               TARGET test_1pni_simple_tpfa
               LABELS porousmediumflow 1p 1pni
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1pni_simple_linear_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pni_simple_linear_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_simple_tpfa params.input -Problem.Name 1pni_simple_linear_tpfa")

dumux_add_test(NAME test_1pni_simple_linear_box
               TARGET test_1pni_simple_box
               LABELS porousmediumflow 1p 1pni
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1pni_simple_linear_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pni_simple_linear_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_simple_box params.input -Problem.Name 1pni_simple_linear_box")

dumux_add_test(NAME test_1pni_simple_source_tpfa
               TARGET test_1pni_simple_tpfa
               LABELS porousmediumflow 1p 1pni
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1pni_simple_source_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pni_simple_source_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_simple_tpfa params.input -Problem.Name 1pni_simple_source_tpfa")

dumux_add_test(NAME test_1pni_simple_source_box
               TARGET test_1pni_simple_box
               LABELS porousmediumflow 1p 1pni
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1pni_simple_source_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pni_simple_source_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_simple_box params.input -Problem.Name 1pni_simple_source_box")
