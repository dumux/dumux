dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME test_md_poromechanics_el2p
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el2p_2p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p_twop-00010.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el2p_poroelastic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p_poroelastic-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p params.input
                       -Problem.Name test_md_poromechanics_el2p")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem_2p.hh
spatialparams_2p.hh
problem_poroelastic.hh
spatialparams_poroelastic.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/poromechanics/el2p)
