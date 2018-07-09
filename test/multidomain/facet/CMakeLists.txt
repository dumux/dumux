dune_symlink_to_source_files(FILES "grid.msh" "test_gridcreator.input")

dune_add_test(NAME test_facetgridcreator_alu
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              SOURCES test_gridcreator.cc
              COMMAND ./test_facetgridcreator_alu
              CMD_ARGS test_gridcreator.input)

dune_add_test(NAME test_facetgridcreator_ug
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-uggrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<3>
              SOURCES test_gridcreator.cc
              COMMAND ./test_facetgridcreator_ug
              CMD_ARGS test_gridcreator.input)

#install sources
install(FILES
test_gridcreator.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/facet)
