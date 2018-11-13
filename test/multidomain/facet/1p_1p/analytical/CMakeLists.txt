dune_symlink_to_source_files(FILES "grids" "params.input" "convergencetest.py")

dune_add_test(NAME test_md_facet_1p1p_tpfa_convergence
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkTpfa LOWDIMTYPETAG=OnePLowDimTpfa
              COMMAND ./convergencetest.py
              CMD_ARGS test_md_facet_1p1p_tpfa_convergence 1e-4)

dune_add_test(NAME test_md_facet_1p1p_box_convergence
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND gmsh_FOUND AND HAVE_UMFPACK )"
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkBox LOWDIMTYPETAG=OnePLowDimBox
              COMMAND ./convergencetest.py
              CMD_ARGS test_md_facet_1p1p_box_convergence 1e-4 )

#install sources
install(FILES
main.cc
problem_bulk.hh
problem_lowdim.hh
spatialparams.hh
convergencetest.py
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet/1p_1p/analytical)
