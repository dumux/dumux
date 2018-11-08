add_input_file_links()

dune_add_test(NAME test_1pncmin_ni_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThermoChemBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pncmin-00064_batch.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pncmin_ni_box-00065.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncmin_ni_box params.input -Problem.Name test_1pncmin_ni_box"
                        --zeroThreshold {"precipVolFrac_CaOH2":5e-6})

#install sources
install(FILES
modifiedcao.hh
modifiedsteamn2cao2h2.hh
problem.hh
spatialparams.hh
reaction.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1pncmin)
