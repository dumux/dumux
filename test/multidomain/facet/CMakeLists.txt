add_subdirectory(1p_1p)
add_subdirectory(1pnc_1pnc)
add_subdirectory(tracer_tracer)

dune_symlink_to_source_files(FILES "grid.msh" "grid2.msh" "test_gridmanager.input" "test_vertexmapper.input" "test_couplingmapper_boundary.input" "2d_grid.msh" "3d_grid.msh")

dumux_add_test(NAME test_facetgridmanager_alu
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              SOURCES test_gridmanager.cc
              COMMAND ./test_facetgridmanager_alu
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetgridmanager_ug
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-uggrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<3>
              SOURCES test_gridmanager.cc
              COMMAND ./test_facetgridmanager_ug
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetcouplingmapper_tpfa_alu
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              SOURCES test_facetcouplingmapper.cc
              COMMAND ./test_facetcouplingmapper_tpfa_alu
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetcouplingmapper_tpfa_ug
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-uggrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<3>
              SOURCES test_facetcouplingmapper.cc
              COMMAND ./test_facetcouplingmapper_tpfa_ug
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetcouplingmapper_mpfa_alu
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              COMPILE_DEFINITIONS USEMPFAINBULK=true
              SOURCES test_facetcouplingmapper.cc
              COMMAND ./test_facetcouplingmapper_mpfa_alu
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetcouplingmapper_mpfa_ug
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-uggrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<3>
              COMPILE_DEFINITIONS USEMPFAINBULK=true
              SOURCES test_facetcouplingmapper.cc
              COMMAND ./test_facetcouplingmapper_mpfa_ug
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetcouplingmapper_box_alu
              LABELS multidomain facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              COMPILE_DEFINITIONS USEBOXINBULK=true
              SOURCES test_facetcouplingmapper.cc
              COMMAND ./test_facetcouplingmapper_box_alu
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetcouplingmapper_box_ug
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-uggrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<3>
              COMPILE_DEFINITIONS USEBOXINBULK=true
              SOURCES test_facetcouplingmapper.cc
              COMMAND ./test_facetcouplingmapper_box_ug
              CMD_ARGS test_gridmanager.input)

dumux_add_test(NAME test_facetcouplingmapper_tpfa_boundary_alu
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
              SOURCES test_facetcouplingmapper_boundary.cc
              COMMAND ./test_facetcouplingmapper_tpfa_boundary_alu
              CMD_ARGS test_couplingmapper_boundary.input)

dumux_add_test(NAME test_facetcouplingmapper_tpfa_boundary_ug
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-uggrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::UGGrid<2>
              SOURCES test_facetcouplingmapper_boundary.cc
              COMMAND ./test_facetcouplingmapper_tpfa_boundary_ug
              CMD_ARGS test_couplingmapper_boundary.input)

dumux_add_test(NAME test_vertexmapper_alu_3d
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>
              COMPILE_DEFINITIONS FACETGRIDTYPE=Dune::FoamGrid<2,3>
              SOURCES test_vertexmapper.cc
              COMMAND ./test_vertexmapper_alu_3d
              CMD_ARGS test_vertexmapper.input -Grid.File 3d_grid.msh)

dumux_add_test(NAME test_vertexmapper_alu_2d
              LABELS multidomain multidomain_facet
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              COMPILE_DEFINITIONS BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
              COMPILE_DEFINITIONS FACETGRIDTYPE=Dune::FoamGrid<1,2>
              SOURCES test_vertexmapper.cc
              COMMAND ./test_vertexmapper_alu_2d
              CMD_ARGS test_vertexmapper.input)
