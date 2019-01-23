add_input_file_links()

dune_add_test(NAME test_1pncminni_discharge_nonequilibrium_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=Discharge_nonequilibrium_Box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pncminni_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1pncminni_box-00065.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pncminni_box params.input -Problem.Name test_1pncminni_box"
                        --zeroThreshold {"precipVolFrac_CaOH2":5e-6})
set(CMAKE_BUILD_TYPE Debug)
#install sources
install(FILES
modifiedcao.hh
problem.hh
spatialparams.hh
reaction.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1pncmin/implicit/discharge_nonequilibrium)
