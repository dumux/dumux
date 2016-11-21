add_input_file_links()

add_dumux_test(test_freeflowstaggered test_freeflowstaggered test_freeflowstaggered.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/1ptestcc-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_staggered")


#install sources
install(FILES
staggeredproblem.hh
test_freeflowstaggered.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
