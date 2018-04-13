dune_symlink_to_source_files(FILES "2d_grid.msh" "3d_grid.msh")

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
