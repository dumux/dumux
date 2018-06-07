add_subdirectory("1p_1p")
dune_symlink_to_source_files(FILES "grid.msh" "grid2.msh" "test_gridcreator.input" "2d_grid.msh" "3d_grid.msh")

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

dune_add_test(NAME test_vertexmapper_alu_3d
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              COMPILE_DEFINITIONS FACETGRIDTYPE=Dune::FoamGrid<2,3>
              SOURCES test_vertexmapper.cc
              COMMAND ./test_vertexmapper_alu_3d
              CMD_ARGS 3d_grid.msh)

dune_add_test(NAME test_vertexmapper_alu_2d
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
              COMPILE_DEFINITIONS FACETGRIDTYPE=Dune::FoamGrid<1,2>
              SOURCES test_vertexmapper.cc
              COMMAND ./test_vertexmapper_alu_2d
              CMD_ARGS 2d_grid.msh)

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
test_gridcreator.cc
test_vertexmapper.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/mixeddimension/facet)
