add_input_file_links(FILES params.input)

dumux_add_test(NAME test_3pni_convection_box
              LABELS porousmediumflow 3p 3pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConvectionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pni_convection_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_3pni_convection_box-00009.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pni_convection_box params.input -Problem.Name test_3pni_convection_box"
                            --zeroThreshold {"velocity_aq \(m/s\)_1":1e-8})

dumux_add_test(NAME test_3pni_convection_tpfa
              LABELS porousmediumflow 3p 3pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConvectionCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pni_convection_tpfa-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_3pni_convection_tpfa-00009.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pni_convection_tpfa params.input -Problem.Name test_3pni_convection_tpfa"
                            --zeroThreshold {"velocity_aq \(m/s\)_1":1e-8})
