dune_symlink_to_source_files(FILES "grids" "params.input" "convergencetest.py")

# executable for tpfa tests
add_executable(test_md_facet_1p1p_tpfa_convergence EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_tpfa_convergence PUBLIC BULKTYPETAG=OnePBulkTpfa LOWDIMTYPETAG=OnePLowDimTpfa)

# exectuable for box tests
add_executable(test_md_facet_1p1p_box_convergence EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_box_convergence PUBLIC BULKTYPETAG=OnePBulkBox LOWDIMTYPETAG=OnePLowDimBox)

dune_add_test(NAME test_md_facet_1p1p_tpfa_convergence
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              TARGET test_md_facet_1p1p_tpfa_convergence
              COMMAND ./convergencetest.py
              CMD_ARGS test_md_facet_1p1p_tpfa_convergence 1e-4)

dune_add_test(NAME test_md_facet_1p1p_box_convergence
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              TARGET test_md_facet_1p1p_box_convergence
              COMMAND ./convergencetest.py
              CMD_ARGS test_md_facet_1p1p_box_convergence 1e-4 )

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
main.cc
problem_bulk.hh
problem_lowdim.hh
spatialparams.hh
convergencetest.py
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1p_1p/analytical)
