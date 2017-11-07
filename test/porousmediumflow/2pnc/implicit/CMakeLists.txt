dune_symlink_to_source_files(FILES test_2pnc.input)

# isothermal tests
dune_add_test(NAME test_box2pnc
              SOURCES test_box2pnc.cc
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pncbox-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/fuelcell-box-00015.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pnc -ParameterFile test_2pnc.input -Problem.Name fuelcell-box")

dune_add_test(NAME test_cc2pnc
              SOURCES test_cc2pnc.cc
              COMMAND python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                             --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pnccc-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/fuelcell-cc-00015.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pnc -ParameterFile test_2pnc.input -Problem.Name fuelcell-cc")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
fuelcellproblem.hh
fuelcellspatialparams.hh
test_box2pnc.cc
test_cc2pnc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pnc)
