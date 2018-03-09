add_input_file_links()

# isothermal tests
add_dumux_test(test_box3p3c test_box3p3c test_box3p3c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/infiltrationbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box3p3c-00005.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box3p3c -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p3c_reference.input")

add_dumux_test(test_cc3p3c test_cc3p3c test_cc3p3c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/infiltrationcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc3p3c-00005.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3c -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p3c_reference.input")

# non-isothermal tests
add_dumux_test(test_box3p3cnicolumnxylol test_box3p3cnicolumnxylol test_box3p3cnicolumnxylol.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/column3p3cnibox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/columnxylolbox-00061.vtu
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
