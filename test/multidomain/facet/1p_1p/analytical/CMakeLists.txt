dune_symlink_to_source_files(FILES "grids" "facetcoupling_1p1p.input" "convergencetest.py")

# executable for tpfa tests
add_executable(test_facetcoupling_tpfa_1p1p EXCLUDE_FROM_ALL test_facetcoupling_fv_1p1p.cc)
target_compile_definitions(test_facetcoupling_tpfa_1p1p PUBLIC BULKTYPETAG=OnePBulkTpfa LOWDIMTYPETAG=OnePLowDimTpfa)

# exectuable for box tests
add_executable(test_facetcoupling_box_1p1p EXCLUDE_FROM_ALL test_facetcoupling_fv_1p1p.cc)
target_compile_definitions(test_facetcoupling_box_1p1p PUBLIC BULKTYPETAG=OnePBulkBox LOWDIMTYPETAG=OnePLowDimBox)

dune_add_test(NAME test_facet_1p1p_tpfa_convergence
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              TARGET test_facetcoupling_tpfa_1p1p
              COMMAND ./convergencetest.py
              CMD_ARGS test_facetcoupling_tpfa_1p1p 1e-4)

dune_add_test(NAME test_facet_1p1p_box_convergence
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              TARGET test_facetcoupling_box_1p1p
              COMMAND ./convergencetest.py
              CMD_ARGS test_facetcoupling_box_1p1p 1e-4)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_facetcoupling_fv_1p1p.cc
bulkProblem.hh
lowdimproblem.hh
spatialparams.hh
convergencetest.py
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1p_1p/analytical)
