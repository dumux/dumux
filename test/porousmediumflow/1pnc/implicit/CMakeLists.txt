add_input_file_links()

dune_add_test(NAME test_1p2c_box
              SOURCES test_1p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/outflowbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnctestbox-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_box test_1p2c_fv.input -Problem.Name 1pnctestbox"
                        --zeroThreshold {"velocity_w \(m/s\)_1":1e-15})

dune_add_test(NAME test_1p2c_tpfa
              SOURCES test_1p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2ctestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnctest_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_tpfa test_1p2c_fv.input -Problem.Name 1pnctest_tpfa"
                        --zeroThreshold {"velocity_w \(m/s\)_1":1e-15})

dune_add_test(NAME test_1p2c_tpfa_restart
              TARGET test_1p2c_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2ctestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnctest_tpfa_restart-00003.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_tpfa test_1p2c_fv.input -Problem.Name 1pnctest_tpfa_restart -TimeLoop.DtInitial 21.4335 -Restart.Time 30.6502 -Restart.File 1pnctest_tpfa-00006.vtu"
                        --zeroThreshold {"velocity_w \(m/s\)_1":1e-15})

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_1p2c_tpfa_restart PROPERTIES DEPENDS test_1p2c_tpfa)

#Note: The mpfa model does not support velocity output yet. We therefore deactivate it for the test and use a copy of the tpfa reference
#      solution where the velocities have been removed.
dune_add_test(NAME test_1p2c_mpfa
              SOURCES test_1p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCTestCCMpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2ctestmpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnctest_mpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_mpfa test_1p2c_fv.input -Problem.Name 1pnctest_mpfa -Vtk.AddVelocity false")

# salt-water intrusion test using the box model
dune_add_test(NAME test_saltwaterintrusion_box
              SOURCES test_1p2c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=SaltWaterIntrusionTestTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/saltwaterintrusion-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/saltwaterintrusion-00035.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_saltwaterintrusion_box saltwaterintrusion.input")

# non-isothermal tests
# conduction
dune_add_test(NAME test_1p2cni_conduction_box
              SOURCES test_1p2cni_conduction_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniboxconduction-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p2cni_conductiontest_box-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_box test_1p2cni_conduction_fv.input -Problem.Name 1p2cni_conductiontest_box"
                        --zeroThreshold {"velocity_w \(m/s\)":1e-9})

dune_add_test(NAME test_1p2cni_conduction_tpfa
              SOURCES test_1p2cni_conduction_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccconduction-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p2cni_conductiontest_tpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_tpfa test_1p2cni_conduction_fv.input -Problem.Name 1p2cni_conductiontest_tpfa"
                        --zeroThreshold {"velocity_w \(m/s\)_0":1e-9})

#Note: The mpfa model does not support velocity output yet. We therefore deactivate it for the test and use a copy of the tpfa reference
#      solution where the velocities have been removed.
dune_add_test(NAME test_1p2cni_conduction_mpfa
              SOURCES test_1p2cni_conduction_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionCCMpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccmpfaconduction-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p2cni_conductiontest_mpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_mpfa test_1p2cni_conduction_fv.input -Problem.Name 1p2cni_conductiontest_mpfa -Vtk.AddVelocity false")

# convection
dune_add_test(NAME test_1p2cni_convection_box
              SOURCES test_1p2cni_convection_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConvectionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniboxconvection-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p2cni_convectiontest_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_box test_1p2cni_convection_fv.input -Problem.Name 1p2cni_convectiontest_box"
                        --zeroThreshold {"velocity_w \(m/s\)":1e-9})

dune_add_test(NAME test_1p2cni_convection_tpfa
              SOURCES test_1p2cni_convection_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConvectionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccconvection-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p2cni_convectiontest_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_tpfa test_1p2cni_convection_fv.input -Problem.Name 1p2cni_convectiontest_tpfa"
                        --zeroThreshold {"velocity_w \(m/s\)_0":1e-9})

#Note: The mpfa model does not support velocity output yet. We therefore deactivate it for the test and use a copy of the tpfa reference
#      solution where the velocities have been removed.
dune_add_test(NAME test_1p2cni_convection_mpfa
              SOURCES test_1p2cni_convection_fv.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConvectionCCMpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccmpfaconvection-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p2cni_convectiontest_mpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_convection_mpfa test_1p2cni_convection_fv.input -Problem.Name 1p2cni_convectiontest_mpfa -Vtk.AddVelocity false")
