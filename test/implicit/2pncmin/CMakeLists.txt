add_input_file_links()

# isothermal tests
add_dumux_test(test_box2pncmin test_box2pncmin test_box2pncmin.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/saltflushbox2pncmin-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/saltflushbox-00044.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pncmin")