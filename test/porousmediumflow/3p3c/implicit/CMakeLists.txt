add_input_file_links()

# isothermal tests
dune_add_test(NAME test_box3p3c
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePThreeCBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/infiltrationbox-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/infiltration3p3cbox-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3p3c test_3p3c_fv.input -Problem.Name infiltration3p3cbox")

dune_add_test(NAME test_cc3p3c
              SOURCES test_3p3c_fv.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePThreeCCCTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/infiltrationcc-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/infiltration3p3ccc-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3c test_3p3c_fv.input -Problem.Name infiltration3p3ccc")

# non-isothermal tests
add_dumux_test(test_box3p3cnicolumnxylol test_box3p3cnicolumnxylol test_box3p3cnicolumnxylol.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/column3p3cnibox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/columnxylolbox-00062.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3p3cnicolumnxylol")

add_dumux_test(test_cc3p3cnicolumnxylol test_cc3p3cnicolumnxylol test_cc3p3cnicolumnxylol.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/column3p3cnicc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/columnxylolcc-00053.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3cnicolumnxylol")

add_dumux_test(test_box3p3cnikuevette test_box3p3cnikuevette test_box3p3cnikuevette.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/kuevette3p3cnibox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/kuevettebox-00004.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3p3cnikuevette")

add_dumux_test(test_cc3p3cnikuevette test_cc3p3cnikuevette test_cc3p3cnikuevette.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/kuevette3p3cnicc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/kuevettecc-00004.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3cnikuevette")

#install sources
install(FILES
columnxylolproblem.hh
columnxylolspatialparams.hh
infiltrationproblem.hh
infiltrationspatialparameters.hh
kuevetteproblem.hh
kuevettespatialparams.hh
test_box3p3c.cc
test_box3p3cnicolumnxylol.cc
test_box3p3cnikuevette.cc
test_cc3p3c.cc
test_cc3p3cnicolumnxylol.cc
test_cc3p3cnikuevette.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/3p3c)
