dune_symlink_to_source_files(FILES "grids" "params.input")

dune_add_test(NAME test_md_facet_1p1p_threedomain_tpfa
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_tpfa_bulk-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_tpfa_facet-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa_facet-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_facet_1p1p_threedomain_tpfa_edge-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa_edge-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_threedomain_tpfa params.input
                       -Problem.Name test_md_facet_1p1p_threedomain_tpfa")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
main.cc
problem_bulk.hh
problem_facet.hh
problem_edge.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1p_1p/threedomain)
