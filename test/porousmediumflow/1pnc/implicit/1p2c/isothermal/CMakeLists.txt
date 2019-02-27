add_subdirectory(saltwaterintrusion)

add_input_file_links()

dune_add_test(NAME test_1p2c_transport_box
              LABELS porousmediumflow 1p2c box
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_transport_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_box params.input -Problem.Name test_1p2c_transport_box"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-15})

dune_add_test(NAME test_1p2c_transport_tpfa
              LABELS porousmediumflow 1p2c tpfa
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_transport_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_tpfa params.input -Problem.Name test_1p2c_transport_tpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-15})

dune_add_test(NAME test_1p2c_transport_mpfa
              LABELS porousmediumflow 1p2c mpfa
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestCCMpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_transport_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_mpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_mpfa params.input -Problem.Name test_1p2c_transport_mpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-13})
