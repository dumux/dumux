add_input_file_links()

# non-isothermal tests
# conduction
dune_add_test(NAME test_1p2cni_conduction_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_conduction_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_box-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_box params.input -Problem.Name test_1p2cni_conduction_box"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-9})

dune_add_test(NAME test_1p2cni_conduction_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_conduction_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_tpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_tpfa params.input -Problem.Name test_1p2cni_conduction_tpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)_0":1e-9})

#Note: The mpfa model does not support velocity output yet. We therefore deactivate it for the test and use a copy of the tpfa reference
#      solution where the velocities have been removed.
dune_add_test(NAME test_1p2cni_conduction_mpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionCCMpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_conduction_mpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_mpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_mpfa params.input -Problem.Name test_1p2cni_conduction_mpfa -Vtk.AddVelocity false")
