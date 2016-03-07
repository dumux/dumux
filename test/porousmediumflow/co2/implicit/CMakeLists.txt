add_input_file_links()

# build target for the CO2 test problem
# Ignore the porosity for all box models since it is defined element-wise in these test
# but the default 2p2c implementation outputs porosity per vertex.
# Depending on the order of the elements, the porosity would differ in these cases.

add_dumux_test(test_boxco2 test_boxco2 test_boxco2.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2box-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousbox-00020.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxco2"
                 --zeroThreshold {"porosity":1})

add_dumux_test(test_ccco2 test_ccco2 test_ccco2.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2cc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneouscc-00019.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccco2")

add_dumux_test(test_restartco2 test_ccco2 test_ccco2.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2box-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousbox-00020.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxco2 -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_restartco2.input"
                 --zeroThreshold {"porosity":1})

# the restart test has to run after the test that produces the restart file
set_tests_properties(test_restartco2 PROPERTIES DEPENDS test_ccco2)

# build target for the CO2 non-isothermal test problem
add_dumux_test(test_boxco2ni test_boxco2ni test_boxco2ni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2nibox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousboxni-00020.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxco2ni"
                 --zeroThreshold {"porosity":1})

add_dumux_test(test_ccco2ni test_ccco2ni test_ccco2ni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2nicc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousccni-00019.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccco2ni")



#install sources
install(FILES
heterogeneousco2tables.hh
heterogeneousproblem.hh
heterogeneousproblemni.hh
heterogeneousspatialparameters.hh
test_boxco2.cc
test_boxco2ni.cc
test_ccco2.cc
test_ccco2ni.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/co2)
