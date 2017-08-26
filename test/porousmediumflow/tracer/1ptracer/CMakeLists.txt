add_input_file_links()

# tracer tests
add_dumux_test(test_1ptracer test_1ptracer test_cctracer.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/cctracer-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cctracer-00027.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cctracer")

#install sources
install(FILES
1ptracertestproblem.hh
1ptracertestspatialparams.hh
test_cctracer.cc
test_cctracer.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/1ptracer)
