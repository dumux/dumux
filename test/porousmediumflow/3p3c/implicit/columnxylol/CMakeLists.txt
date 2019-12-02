add_input_file_links(FILES params.input)

# non-isothermal tests
dumux_add_test(NAME test_3p3cni_columnxylol_box
              LABELS porousmediumflow 3p3c 3p3cni
              SOURCES main.cc
              TIMEOUT 1500
              COMPILE_DEFINITIONS TYPETAG=ColumnBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3cni_columnxylol_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_columnxylol_box-00069.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_columnxylol_box params.input -Problem.Name test_3p3cni_columnxylol_box")

dumux_add_test(NAME test_3p3cni_columnxylol_tpfa
              LABELS porousmediumflow 3p3c 3p3cni
              SOURCES main.cc
              TIMEOUT 1500
              COMPILE_DEFINITIONS TYPETAG=ColumnCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3cni_columnxylol_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_columnxylol_tpfa-00053.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3cni_columnxylol_tpfa params.input -Problem.Name test_3p3cni_columnxylol_tpfa")
