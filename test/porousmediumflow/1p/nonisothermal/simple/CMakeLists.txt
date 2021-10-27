add_input_file_links()

dumux_add_test(NAME test_1pni_without_divergence_tpfa
              LABELS porousmediumflow 1p 1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNISimpleTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pni_without_divergence_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pni_without_divergence_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_without_divergence_tpfa params.input -Problem.Name 1pni_without_divergence_tpfa")

dumux_add_test(NAME test_1pni_without_divergence_box
              LABELS porousmediumflow 1p 1pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNISimpleBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/1pni_without_divergence_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1pni_without_divergence_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_without_divergence_box params.input -Problem.Name 1pni_without_divergence_box")
