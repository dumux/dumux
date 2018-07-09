dune_symlink_to_source_files(FILES "el1p.input")

dune_add_test(NAME test_el1p
              SOURCES test_el1p.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/el1p_p_reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/onep.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/el1p_u_reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/poromech.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_el1p el1p.input")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
1pproblem.hh
1pspatialparams.hh
poroelasticproblem.hh
poroelasticspatialparams.hh
test_el1p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/poromechanics/el1p)
