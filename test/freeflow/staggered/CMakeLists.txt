add_input_file_links()

add_dumux_test(test_cc1p test_cc1p test_cc1p.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc1p")


#install sources
install(FILES
1ptestproblem.hh
1ptestspatialparams.hh
test_cc1p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
