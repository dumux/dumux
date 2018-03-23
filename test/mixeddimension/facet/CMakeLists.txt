dune_symlink_to_source_files(FILES "grid.msh" "test_gridcreator.input")

dune_add_test(NAME test_facetgridcreator_alu
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              SOURCES test_gridcreator.cc
              COMMAND ./test_facetgridcreator_alu
              CMD_ARGS test_gridcreator.input)

dune_add_test(NAME test_facetgridcreator_ug
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<3>
              SOURCES test_gridcreator.cc
              COMMAND ./test_facetgridcreator_ug
              CMD_ARGS test_gridcreator.input)

dune_add_test(NAME test_tpfafacetcouplingmapper_alu
             CMAKE_GUARD dune-foamgrid_FOUND
             CMAKE_GUARD dune-alugrid_FOUND
             COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
             SOURCES test_tpfafacetcouplingmapper.cc
             COMMAND ./test_tpfafacetcouplingmapper_alu
             CMD_ARGS test_gridcreator.input)

dune_add_test(NAME test_tpfafacetcouplingmapper_ug
             CMAKE_GUARD dune-foamgrid_FOUND
             CMAKE_GUARD dune-uggrid_FOUND
             COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<3>
             SOURCES test_tpfafacetcouplingmapper.cc
             COMMAND ./test_tpfafacetcouplingmapper_ug
             CMD_ARGS test_gridcreator.input)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_gridcreator.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/mixeddimension/facet)
