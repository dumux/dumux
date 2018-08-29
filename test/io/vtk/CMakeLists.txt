dune_add_test(NAME test_vtkreader_3d
              SOURCES test_vtkreader.cc
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::cube,Dune::nonconforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_3d
                                  ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p2c_richards2c_3d-reference.vtu test-3d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1d3d_1p2c_richards2c_3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test-3d.vtu)

# dune_add_test(NAME test_vtkreader_3d_structured_notimplemented
#               SOURCES test_vtkreader.cc
#               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>
#               EXPECT_COMPILE_FAIL)

dune_add_test(NAME test_vtkreader_2d
              SOURCES test_vtkreader.cc
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_2d
                                  ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu test-2d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/1ptestbox-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test-2d.vtu)

dune_add_test(NAME test_vtkreader_2d3d
              SOURCES test_vtkreader.cc
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,3,Dune::simplex,Dune::conforming>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_2d3d
                                  ${CMAKE_SOURCE_DIR}/test/references/test_embedded_fracture_tpfa_2d.vtu test-2d3d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_fracture_tpfa_2d.vtu
                   ${CMAKE_CURRENT_BINARY_DIR}/test-2d3d.vtu)

dune_add_test(NAME test_vtkreader_1d
              SOURCES test_vtkreader.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::FoamGrid<1,3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_vtkreader_1d
                                  ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1p_richards_tpfa_1d.vtp test-1d"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_embedded_1p_richards_tpfa_1d.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test-1d.vtp)

dune_add_test(NAME test_staggered_ff_pv_names_1d
             COMPILE_DEFINITIONS DIM=1
             SOURCES test_staggeredfreeflowpvnames.cc)

dune_add_test(NAME test_staggered_ff_pv_names_2d
             COMPILE_DEFINITIONS DIM=2
             SOURCES test_staggeredfreeflowpvnames.cc)

dune_add_test(NAME test_staggered_ff_pv_names_3d
             COMPILE_DEFINITIONS DIM=3
             SOURCES test_staggeredfreeflowpvnames.cc)
