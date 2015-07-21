# build target for the CO2 test problem
add_dumux_test(test_boxco2 test_boxco2 test_boxco2.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2box-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousbox-00020.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxco2")

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
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxco2 -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_restartco2.input")

  # build target for the CO2 non-isothermal test problem
add_dumux_test(test_boxco2ni test_boxco2ni test_boxco2ni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2nibox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousboxni-00020.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxco2ni")

add_dumux_test(test_ccco2ni test_ccco2ni test_ccco2ni.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/co2nicc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/heterogeneousccni-00019.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccco2ni")


