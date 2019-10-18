add_input_file_links()

# non-isothermal tests
# conduction
dumux_add_test(NAME test_1p2c_nonequilibrium_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCThermalNonequilibriumBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_nonequilibrium_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_nonequilibrium_box-00045.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_nonequilibrium_box params.input -Problem.Name test_1p2c_nonequilibrium_box"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-15})

dumux_add_test(NAME test_1p2c_nonequilibrium_tpfa
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCThermalNonequilibriumCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_nonequilibrium_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_nonequilibrium_tpfa-00044.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_nonequilibrium_tpfa params.input -Problem.Name test_1p2c_nonequilibrium_tpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_1":1e-15})
