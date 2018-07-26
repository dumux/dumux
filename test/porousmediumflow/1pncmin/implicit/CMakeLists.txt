add_subdirectory("nonisothermal")

add_input_file_links()

dune_add_test(NAME test_1pncminni_discharge_box
              SOURCES test_1pncminni_discharge_fv.cc
              COMPILE_DEFINITIONS TYPETAG=DischargeBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pncmin-00064_batch.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pncmintest_batch_box-00065.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncminni_box test_1pncminni_discharge.input -Problem.Name 1pncmintest_batch_box"
                        --zeroThreshold {"precipVolFrac_CaOH2":5e-6})

dune_add_test(NAME test_1pncminni_charge_box
              SOURCES test_1pncminni_charge_fv.cc
              COMPILE_DEFINITIONS TYPETAG=ChargeBoxTypeTag
              CMD_ARGS test_1pncminni_charge.input)

#install sources
install(FILES
modifiedcao.hh
charge_problem.hh
discharge_problem.hh
thermochemspatialparams.hh
thermochemreaction.hh
test_1pncminni_charge_fv.cc
test_1pncminni_discharge_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1pncmin)
>>>>>>> d228c01c1... [1pncmin] split up test into charge and discharge
