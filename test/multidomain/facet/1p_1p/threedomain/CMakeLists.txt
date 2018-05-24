dune_symlink_to_source_files(FILES "grids" "test_facetcoupling_1p1p_threedomain.input")

dune_add_test(NAME test_facetcoupling_tpfa_1p1p_threedomain
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              SOURCES test_facetcoupling_tpfa_1p1p_threedomain.cc
              COMMAND ./test_facetcoupling_tpfa_1p1p_threedomain
              CMD_ARGS test_facetcoupling_1p1p_threedomain.input)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_facetcoupling_tpfa_1p1p_threedomain.cc
bulkProblem.hh
facetproblem.hh
edgeproblem.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1p_1p/threedomain)
