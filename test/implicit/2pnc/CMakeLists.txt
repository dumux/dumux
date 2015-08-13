add_input_file_links()

# isothermal tests
add_dumux_test(test_box2pnc test_box2pnc test_box2pnc.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pnc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/fuelcell-00016.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pnc")