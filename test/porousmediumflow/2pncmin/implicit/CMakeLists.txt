dune_symlink_to_source_files(FILES test_2pncmin.input)

# isothermal tests
dune_add_test(NAME test_2pncmin_box
              SOURCES test_2pncmin_fv.cc
              COMPILE_DEFINITIONS TYPETAG=DissolutionBoxTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/saltflushbox2pncmin-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/saltflushbox-00050.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_box -ParameterFile test_2pncmin.input -Problem.Name saltflushbox")

dune_add_test(NAME test_2pncmin_tpfa
              SOURCES test_2pncmin_fv.cc
              COMPILE_DEFINITIONS TYPETAG=DissolutionCCTpfaTypeTag
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/saltflushcc2pncmin-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/saltflushtpfa-00043.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_tpfa -ParameterFile test_2pncmin.input -Problem.Name saltflushtpfa")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
dissolutionproblem.hh
dissolutionspatialparams.hh
test_2pncmin_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pncmin)
