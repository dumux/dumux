dune_symlink_to_source_files(FILES "poroelastic.input")

# using box and numeric differentiation
dune_add_test(NAME test_poroelastic_box
              SOURCES test_poroelastic.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/poroelasticbox-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/poroelastic-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_poroelastic_box poroelastic.input")

set(CMAKE_BUILD_TYPE Release)

install(FILES
problem.hh
spatialparams.hh
test_poroelastic.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test/geomechanics/poroelastic)
