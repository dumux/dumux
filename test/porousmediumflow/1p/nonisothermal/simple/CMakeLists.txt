add_input_file_links()

dumux_add_test(NAME test_1pni_simple_tpfa
              LABELS porousmediumflow 1p 1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNISimpleTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa params.input -Problem.Name test_1pni_conduction_mpfa")

dumux_add_test(NAME test_1pni_simple_box
              LABELS porousmediumflow 1p 1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNISimpleBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa params.input -Problem.Name test_1pni_conduction_mpfa")
