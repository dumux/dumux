#add links to input files
add_input_file_links()

add_dumux_test(test_boxsteaminjection test_boxsteaminjection test_boxsteaminjection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectionbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_boxsteaminjection-00008.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxsteaminjection")

add_dumux_test(test_ccsteaminjection test_ccsteaminjection test_ccsteaminjection.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/steaminjectioncc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_ccsteaminjection-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccsteaminjection")
