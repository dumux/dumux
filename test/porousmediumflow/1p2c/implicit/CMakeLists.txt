add_input_file_links()

# isothermal tests
add_dumux_test(test_box1p2c test_box1p2c test_box1p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/outflowbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/outflowbox-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2c"
                 --zeroThreshold {"velocity":5e-16})

add_dumux_test(test_cc1p2c test_cc1p2c test_cc1p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2ctestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1p2ctestcc-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2c"
                 --zeroThreshold {"velocity":5e-16})

add_dumux_test(test_ccmpfa1p2c test_ccmpfa1p2c test_ccmpfa1p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2ctestccmpfa-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1p2ctestccmpfa-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpfa1p2c"
                 --zeroThreshold {"velocity":5e-16})

# non-isothermal tests
add_dumux_test(test_box1p2cniconduction test_box1p2cniconduction test_box1p2cniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniboxconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1p2cniconduction-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2cniconduction"
                 --zeroThreshold {"velocity":1e-8})

add_dumux_test(test_box1p2cniconvection test_box1p2cniconvection test_box1p2cniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniboxconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1p2cniconvection-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2cniconvection"
                 --zeroThreshold {"velocity":1e-16})

add_dumux_test(test_cc1p2cniconduction test_cc1p2cniconduction test_cc1p2cniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1p2cniconduction-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2cniconduction"
                 --zeroThreshold {"velocity":1e-8})

add_dumux_test(test_cc1p2cniconvection test_cc1p2cniconvection test_cc1p2cniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1p2cniconvection-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2cniconvection")

add_dumux_test(test_ccmpfa1p2cniconduction test_ccmpfa1p2cniconduction test_ccmpfa1p2cniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccmpfaconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/ccmpfa1p2cniconduction-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpfa1p2cniconduction")

add_dumux_test(test_ccmpfa1p2cniconvection test_ccmpfa1p2cniconvection test_ccmpfa1p2cniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccmpfaconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/ccmpfa1p2cniconvection-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpfa1p2cniconvection")

#install sources
install(FILES
1p2cniconductionproblem.hh
1p2cniconvectionproblem.hh
1p2cnispatialparams.hh
1p2coutflowproblem.hh
1p2coutflowspatialparams.hh
test_box1p2c.cc
test_box1p2cniconduction.cc
test_box1p2cniconvection.cc
test_cc1p2c.cc
test_cc1p2cniconduction.cc
test_cc1p2cniconvection.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p2c)
