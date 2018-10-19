dune_symlink_to_source_files(FILES "grids" "facetcoupling_1pnc1pnc.input")

dune_add_test(NAME test_facet_1pnc1pnc_tpfa
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkTpfa
              COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePNCLowDimTpfa
              SOURCES test_facetcoupling_fv_1pnc_1pnc.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/facetcoupling_1pnc1pnc_bulk_tpfa.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnc_1pnc_bulk_tpfa-00016.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/facetcoupling_1pnc1pnc_lowdim_tpfa.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnc_1pnc_lowdim_tpfa-00016.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_facet_1pnc1pnc_tpfa facetcoupling_1pnc1pnc.input -Bulk.Problem.Name 1pnc_1pnc_bulk_tpfa -LowDim.Problem.Name 1pnc_1pnc_lowdim_tpfa")

dune_add_test(NAME test_facet_1pnc1pnc_box
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              COMPILE_DEFINITIONS BULKTYPETAG=OnePNCBulkBox
              COMPILE_DEFINITIONS LOWDIMTYPETAG=OnePNCLowDimBox
              SOURCES test_facetcoupling_fv_1pnc_1pnc.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/facetcoupling_1pnc1pnc_bulk_box.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnc_1pnc_bulk_box-00016.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/facetcoupling_1pnc1pnc_lowdim_box.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/1pnc_1pnc_lowdim_box-00016.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_facet_1pnc1pnc_box facetcoupling_1pnc1pnc.input -Bulk.Problem.Name 1pnc_1pnc_bulk_box -LowDim.Problem.Name 1pnc_1pnc_lowdim_box")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_facetcoupling_fv_1pnc1pnc.cc
bulkProblem.hh
lowdimproblem.hh
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1pnc_1pnc/)
