add_input_file_links(FILES params.input)

dumux_add_test(NAME test_3p3cni_kuevette_box
              LABELS porousmediumflow 3p3c 3p3cni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=KuevetteBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3cni_kuevette_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_box-00004.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_box params.input -Problem.Name test_3p3cni_kuevette_box")

dumux_add_test(NAME test_3p3cni_kuevette_tpfa
              LABELS porousmediumflow 3p3c 3p3cni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=KuevetteCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3cni_kuevette_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_tpfa-00004.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_kuevette_tpfa params.input -Problem.Name test_3p3cni_kuevette_tpfa")
