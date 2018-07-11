dune_symlink_to_source_files(FILES "grids" "facetcoupling_1p1p_threedomain.input")

dune_add_test(NAME test_facetcoupling_tpfa_1p1p_threedomain
              SOURCES test_facetcoupling_tpfa_1p1p_threedomain.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCCTypeTag LOWDIMTYPETAG=BloodFlowCCTypeTag COUPLINGMODE=EmbeddedCouplingMode::average
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/facetcoupling_1p1p_threedomain_bulk.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_1p_bulk-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/facetcoupling_1p1p_threedomain_facet.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_1p_facet-00001.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/facetcoupling_1p1p_threedomain_edge.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_1p_edge-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_facetcoupling_tpfa_1p1p_threedomain facetcoupling_1p1p_threedomain.input")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_facetcoupling_tpfa_1p1p_threedomain.cc
bulkProblem.hh
facetproblem.hh
edgeproblem.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1p_1p/threedomain)
