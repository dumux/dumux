add_input_file_links()

add_dumux_test(test_stokes test_stokes test_stokes.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes-00013.vtu)

#install sources
install(FILES
stokestestproblem.hh
test_stokes.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/stokes)
