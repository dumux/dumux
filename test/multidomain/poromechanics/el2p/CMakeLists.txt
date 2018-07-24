dune_symlink_to_source_files(FILES "el2p.input")

dune_add_test(NAME test_el2p
              SOURCES test_el2p.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/el2p_2p_reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/twop-00016.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/el2p_u_reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/poroelastic-00016.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_el2p el2p.input")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
2pproblem.hh
2pspatialparams.hh
poroelasticproblem.hh
poroelasticspatialparams.hh
test_el2p.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/poromechanics/el2p)
