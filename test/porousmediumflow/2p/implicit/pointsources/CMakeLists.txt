add_dumux_test(test_ccadaptive2ppointsource test_ccadaptive2ppointsource test_ccadaptive2ppointsource.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/lensccadaptivepointsource-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/lensccadaptive-00012.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccadaptive2ppointsource")

#install sources
install(FILES
lensproblem.hh
test_ccadaptive2ppointsource.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/pointsources)
