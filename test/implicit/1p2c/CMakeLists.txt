# isothermal tests
add_dumux_test(test_box1p2c test_box1p2c test_box1p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/outflowbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/outflowbox-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2c")

add_dumux_test(test_cc1p2c test_cc1p2c test_cc1p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/outflowcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/outflowcc-00010.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2c")

# non-isothermal tests
add_dumux_test(test_box1p2cniconduction test_box1p2cniconduction test_box1p2cniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniboxconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1p2cniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2cniconduction")

add_dumux_test(test_box1p2cniconvection test_box1p2cniconvection test_box1p2cniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniboxconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/box1p2cniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2cniconvection")

add_dumux_test(test_cc1p2cniconduction test_cc1p2cniconduction test_cc1p2cniconduction.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccconduction-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1p2cniconduction-00007.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2cniconduction")

add_dumux_test(test_cc1p2cniconvection test_cc1p2cniconvection test_cc1p2cniconvection.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1p2cniccconvection-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cc1p2cniconvection-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2cniconvection")
