add_input_file_links()

# isothermal tests
add_dumux_test(lensbox lensbox lensbox.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/satcondition-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/satcondition-00009.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test1d_satcondition")

#install sources
install(FILES
lensboxproblem.hh
lensboxparams.hh
lensbox.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/lensbox)
