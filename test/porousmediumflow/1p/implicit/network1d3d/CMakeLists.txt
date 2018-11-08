add_input_file_links()
add_gstat_file_links()

dune_symlink_to_source_files(FILES grids)
dune_symlink_to_source_files(FILES "convergencetest.py")

# dim < dimWorld tests with Dune::Foamgrid<1,3>
dune_add_test(NAME test_1p_network1d3d_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestCCTpfa
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_network1d3d_tpfa params.input -Problem.Name test_1p_network1d3d_tpfa)

dune_add_test(NAME test_1p_network1d3d_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TubesTestBox
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_network1d3d_box params.input -Problem.Name test_1p_network1d3d_box)

#install sources
install(FILES
problem.hh
spatialparams.hh
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/1p)
