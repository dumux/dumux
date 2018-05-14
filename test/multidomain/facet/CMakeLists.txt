add_subdirectory("1p_1p")
add_subdirectory("box")
dune_symlink_to_source_files(FILES "grid.msh" "grid2.msh" "test_gridcreator.input")

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

dune_add_test(NAME test_tpfafacetcouplingmapper_boundary_alu
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
              SOURCES test_tpfafacetcouplingmapper_boundary.cc
              COMMAND ./test_tpfafacetcouplingmapper_boundary_alu)

dune_add_test(NAME test_tpfafacetcouplingmapper_boundary_ug
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-uggrid_FOUND
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<2>
              SOURCES test_tpfafacetcouplingmapper_boundary.cc
              COMMAND ./test_tpfafacetcouplingmapper_boundary_ug)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_gridcreator.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/mixeddimension/facet)
