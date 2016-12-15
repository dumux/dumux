add_input_file_links()

add_dumux_test(test_liddrivencavity test_closedsystem test_closedsystem.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/liddrivencavity-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_liddrivencavity-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_liddrivencavity.input")

add_dumux_test(test_hydrostaticpressure test_closedsystem test_closedsystem.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/stokeshydrostaticpressure-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/test_hydrostaticpressure-00006.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_closedsystem test_hydrostaticpressure.input")


#install sources
install(FILES
staggeredproblem.hh
test_freeflowstaggered.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
