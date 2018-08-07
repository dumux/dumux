add_input_file_links()

# isothermal tests
dune_add_test(NAME test_box3p
              SOURCES test_3p_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/infiltration3pbox-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pbox-00008.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3p test_3p_fv.input -Problem.Name infiltration3pbox")

dune_add_test(NAME test_cc3p
              SOURCES test_3p_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/infiltration3pcc-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pcc-00008.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3p test_3p_fv.input -Problem.Name infiltration3pcc")

# non-isothermal tests
dune_add_test(NAME test_box3pniconvection
              SOURCES test_3pni_fv_convection.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConvectionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/3pniconvectionbox-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconvection-00009.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconvection test_3pni_fv_convection.input -Problem.Name test_box3pniconvection"
                            --zeroThreshold {"velocity_w \(m/s\)_1":1e-8})

dune_add_test(NAME test_cc3pniconvection
              SOURCES test_3pni_fv_convection.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConvectionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/3pniconvectioncc-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconvection-00009.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconvection test_3pni_fv_convection.input -Problem.Name test_cc3pniconvection"
                            --zeroThreshold {"velocity_w \(m/s\)_1":1e-8})

dune_add_test(NAME test_box3pniconduction
              SOURCES test_3pni_fv_conduction.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConductionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/3pniconductionbox-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconduction-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3pniconduction test_3pni_fv_conduction.input -Problem.Name test_box3pniconduction"
                            --zeroThreshold {"velocity_w \(m/s\)_1":1e-8})

dune_add_test(NAME test_cc3pniconduction
              SOURCES test_3pni_fv_conduction.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConductionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/3pniconductioncc-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconduction-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3pniconduction test_3pni_fv_conduction.input -Problem.Name test_cc3pniconduction"
                            --zeroThreshold {"velocity_w \(m/s\)_1":1e-8})

#install sources
install(FILES
3pniconductionproblem.hh
3pniconvectionproblem.hh
3pnispatialparams.hh
infiltration3pproblem.hh
infiltration3pspatialparams.hh
test_box3p.cc
test_box3pniconduction.cc
test_box3pniconvection.cc
test_cc3p.cc
test_cc3pniconduction.cc
test_cc3pniconvection.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3p)
