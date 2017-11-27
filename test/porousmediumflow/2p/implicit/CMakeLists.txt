#add_subdirectory(pointsources)
add_subdirectory(fracture)
add_subdirectory(incompressible)
add_subdirectory(nonisothermal)
add_input_file_links()
dune_symlink_to_source_files(FILES grids)

# # isothermal tests
# add_dumux_test(test_ccadaptive2p test_ccadaptive2p test_ccadaptive2p.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/lensccadaptive-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/lensccadaptive-00017.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_ccadaptive2p")
#
# add_dumux_test(test_boxadaptive2p test_boxadaptive2p test_boxadaptive2p.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/lensboxadaptive-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/lensboxadaptive-00014.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_boxadaptive2p")
#
# if(HAVE_OPM_GRID)
# add_dumux_test(test_cc2pcornerpoint test_cc2pcornerpoint test_cc2pcornerpoint.cc
#                python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
#                  --script fuzzy
#                  --files ${CMAKE_SOURCE_DIR}/test/references/cc2pcornerpoint-reference.vtu
#                          ${CMAKE_CURRENT_BINARY_DIR}/cc2pcornerpoint-00005.vtu
#                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pcornerpoint")
# endif()

#install sources
install(FILES
cc2pcornerpointproblem.hh
cc2pcornerpointspatialparams.hh
generalizeddirichletproblem.hh
generalizeddirichletspatialparams.hh
test_boxadaptive2p.cc
test_cc2pcornerpoint.cc
test_ccadaptive2p.cc
test_generalizeddirichlet.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p)
