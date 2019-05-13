add_input_file_links()
dune_symlink_to_source_files(FILES polyline.vtp)

dune_add_test(NAME test_vtkreader_3d
              SOURCES test_vtkreader.cc
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::cube,Dune::nonconforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_3d
                                  ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p2c_richards2c_3d-reference.vtu test-3d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p2c_richards2c_3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test-3d.vtu)

# dune_add_test(NAME test_vtkreader_3d_structured_notimplemented
#               SOURCES test_vtkreader.cc
#               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>
#               EXPECT_COMPILE_FAIL)

dune_add_test(NAME test_vtkreader_2d
              SOURCES test_vtkreader.cc
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_2d
                                  ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu test-2d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test-2d.vtu)

dune_add_test(NAME test_vtkreader_2d3d
              SOURCES test_vtkreader.cc
              LABELS unit io
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::FoamGrid<2,3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_2d3d
                                  ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded2d3d_fracture1p_tpfa_2d-reference.vtu test-2d3d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded2d3d_fracture1p_tpfa_2d-reference.vtu
                   ${CMAKE_CURRENT_BINARY_DIR}/test-2d3d.vtu)

dune_add_test(NAME test_vtkreader_1d
              SOURCES test_vtkreader.cc
              LABELS unit io
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::FoamGrid<1,3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_1d
                                  ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp test-1d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded1d3d_1p_richards_tpfatpfa_1d-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test-1d.vtp)

dune_add_test(NAME test_vtkreader_1d_polyline
              TARGET test_vtkreader_1d
              LABELS unit io
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_1d polyline.vtp test_polyline"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_vtkreader_1d_polyline-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test_polyline.vtp)

dune_add_test(NAME test_vtk_staggeredfreeflowpvnames
              SOURCES test_vtk_staggeredfreeflowpvnames.cc
              LABELS unit io)
