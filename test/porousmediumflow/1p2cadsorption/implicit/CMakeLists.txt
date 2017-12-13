add_input_file_links()

# isothermal tests
add_dumux_test(test_box1p2cadsorption test_box1p2cadsorption test_box1p2cadsorption.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/outflowbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/outflowbox-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p2cadsorption"
                 --zeroThreshold {"velocity":5e-16})

add_dumux_test(test_cc1p2cadsorption test_cc1p2cadsorption test_cc1p2cadsorption.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/outflowcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/outflowcc-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2cadsorption"
                 --zeroThreshold {"velocity":5e-16})

#install sources
install(FILES
1p2cadsorptionoutflowproblem.hh
1p2cadsorptionoutflowspatialparams.hh
test_box1p2cadsorption.cc
test_cc1p2cadsorption.cc
test_1p2cadsorption.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p2cadsorption)
