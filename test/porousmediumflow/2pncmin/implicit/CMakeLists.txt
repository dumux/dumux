# isothermal tests
add_dumux_test(test_box2pncmin test_box2pncmin test_box2pncmin.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/saltflushbox2pncmin-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/saltflushbox-00048.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pncmin -ParameterFile test_2pncmin.input -Problem.Name saltflushbox")

add_dumux_test(test_cc2pncmin test_cc2pncmin test_cc2pncmin.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/saltflushcc2pncmin-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/saltflushcc-00043.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pncmin -ParameterFile test_2pncmin.input -Problem.Name saltflushcc")

dune_symlink_to_source_files(FILES "test_2pncmin.input")

#install sources
install(FILES
dissolutionproblem.hh
dissolutionspatialparams.hh
test_box2pncmin.cc
test_cc2pncmin.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pncmin)
