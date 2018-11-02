dune_symlink_to_source_files(FILES test_2pncmin.input)

# isothermal tests
dune_add_test(NAME test_2pncmin_box
              SOURCES test_2pncmin_fv.cc
              COMPILE_DEFINITIONS TYPETAG=DissolutionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/saltflushbox2pncmin-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/saltflushbox-00044.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_box -ParameterFile test_2pncmin.input -Problem.Name saltflushbox")

dune_add_test(NAME test_2pncmin_box_restart
              TARGET  test_2pncmin_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/saltflushbox2pncmin-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/saltflushbox_restart-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_box test_2pncmin.input -Problem.Name saltflushbox_restart -TimeLoop.DtInitial 50000 -Restart.Time 756290 -Restart.File saltflushbox-00039.vtu")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_2pncmin_box_restart PROPERTIES DEPENDS test_2pncmin_box)

dune_add_test(NAME test_2pncmin_tpfa
              SOURCES test_2pncmin_fv.cc
              COMPILE_DEFINITIONS TYPETAG=DissolutionCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/saltflushtpfa2pncmin-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/saltflushtpfa-00043.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_tpfa -ParameterFile test_2pncmin.input -Problem.Name saltflushtpfa")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
dissolutionproblem.hh
dissolutionspatialparams.hh
test_2pncmin_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pncmin)
