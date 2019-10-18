add_input_file_links()

dumux_add_test(NAME test_1pncminni_box
              LABELS porousmediumflow 1pncmin  1pncminni
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=ThermoChemBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pncminni_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pncminni_box-00065.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncminni_box params.input -Problem.Name test_1pncminni_box"
                        --zeroThreshold {"precipVolFrac_CaOH2":5e-6})
