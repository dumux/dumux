# isothermal tests
add_dumux_test(test_box2pnc test_box2pnc test_box2pnc.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pncbox-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/fuelcell-box-00016.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_box2pnc -ParameterFile test_2pnc.input -Problem.Name fuelcell-box")

add_dumux_test(test_cc2pnc test_cc2pnc test_cc2pnc.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/fuelcell2pnccc-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/fuelcell-cc-00016.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_cc2pnc -ParameterFile test_2pnc.input -Problem.Name fuelcell-cc")

dune_symlink_to_source_files(FILES test_2pnc.input)

#install sources
install(FILES
fuelcellproblem.hh
fuelcellspatialparams.hh
test_box2pnc.cc
test_cc2pnc.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pnc)
