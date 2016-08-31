#add links to input files
add_input_file_links()

add_dumux_test(test_boxsteaminjection test_boxsteaminjection test_boxsteaminjection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestbox-00001.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box1p")
