add_input_file_links()

# tracer tests
add_dumux_test(test_cctracer test_cctracer test_cctracer.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/cctracer-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/cctracer-00027.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cctracer")

add_dumux_test(test_ccmpfatracer test_ccmpfatracer test_ccmpfatracer.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/ccmpfatracer-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/ccmpfatracer-00027.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccmpfatracer")

add_dumux_test(test_boxtracer test_boxtracer test_boxtracer.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/boxtracer-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/boxtracer-00027.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxtracer")

#install sources
install(FILES
tracertestproblem.hh
tracertestspatialparams.hh
test_cctracer.cc
test_cctracer.input
test_ccmpfatracer.cc
test_ccmpfatracer.input
test_boxtracer.cc
test_boxtracer.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/implicit)
