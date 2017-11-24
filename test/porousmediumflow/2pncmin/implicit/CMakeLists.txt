dune_symlink_to_source_files(FILES test_2pncmin.input)

# isothermal tests
dune_add_test(NAME test_box2pncmin
              SOURCES test_box2pncmin.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/saltflushbox2pncmin-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/saltflushbox-00046.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pncmin -ParameterFile test_2pncmin.input -Problem.Name saltflushbox")

dune_add_test(NAME test_cc2pncmin
              SOURCES test_cc2pncmin.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/saltflushcc2pncmin-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/saltflushcc-00065.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pncmin -ParameterFile test_2pncmin.input -Problem.Name saltflushcc")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
dissolutionproblem.hh
dissolutionspatialparams.hh
test_box2pncmin.cc
test_cc2pncmin.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pncmin)
