add_input_file_links()

# non-isothermal tests
# convection
dumux_add_test(NAME test_1p2cni_convection_box
              LABELS porousmediumflow 1p2c  1p2cni
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConvectionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_convection_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_box params.input -Problem.Name test_1p2cni_convection_box"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-9})

dumux_add_test(NAME test_1p2cni_convection_tpfa
              LABELS porousmediumflow 1p2c  1p2cni
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConvectionCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_convection_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_tpfa params.input -Problem.Name test_1p2cni_convection_tpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_0":1e-9})

dumux_add_test(NAME test_1p2cni_convection_mpfa
              LABELS porousmediumflow 1p2c  1p2cni
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConvectionCCMpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_convection_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_mpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_mpfa params.input -Problem.Name test_1p2cni_convection_mpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_0":1e-9})
