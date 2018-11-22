add_subdirectory(saltwaterintrusion)

add_input_file_links()

dune_add_test(NAME test_1p2c_transport_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_transport_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_box params.input -Problem.Name test_1p2c_transport_box"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-15})

dune_add_test(NAME test_1p2c_transport_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_transport_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_tpfa params.input -Problem.Name test_1p2c_transport_tpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-15})

#Note: The mpfa model does not support velocity output yet. We therefore deactivate it for the test and use a copy of the tpfa reference
#      solution where the velocities have been removed.
dune_add_test(NAME test_1p2c_transport_mpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestCCMpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_transport_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_mpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_transport_mpfa params.input -Problem.Name test_1p2c_transport_mpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-13})
