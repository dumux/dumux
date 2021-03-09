dune_symlink_to_source_files(FILES "params.input" "grids")

# quadrilaterals alu grid
dumux_add_test(NAME test_2p_boxdfm_quads_alu
              LABELS porousmediumflow 2p boxdfm
              SOURCES main.cc
              CMAKE_GUARD dune-alugrid_FOUND
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_boxdfm_2d_quads-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_quads_alu-00023.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_quads_alu params.input -Problem.Name test_2p_boxdfm_quads_alu -Grid.File grids/quads_2d.msh")

# quadrilaterals ug grid
dumux_add_test(NAME test_2p_boxdfm_quads_ug
              LABELS porousmediumflow 2p boxdfm
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_boxdfm_2d_quads-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_quads_ug-00023.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_quads_ug params.input -Problem.Name test_2p_boxdfm_quads_ug -Grid.File grids/quads_2d.msh")

# triangles alu
dumux_add_test(NAME test_2p_boxdfm_trias_alu
              LABELS porousmediumflow 2p boxdfm
              TIMEOUT 1500
             SOURCES main.cc
             CMAKE_GUARD dune-foamgrid_FOUND
             CMAKE_GUARD dune-alugrid_FOUND
             COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
             COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_boxdfm_2d_trias-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_trias_alu-00064.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_trias_alu params.input -Problem.Name test_2p_boxdfm_trias_alu -Grid.File grids/durlofsky.msh")

# triangles ug
dumux_add_test(NAME test_2p_boxdfm_trias_ug
              LABELS porousmediumflow 2p boxdfm
              TIMEOUT 1500
             SOURCES main.cc
             CMAKE_GUARD dune-uggrid_FOUND
             CMAKE_GUARD dune-foamgrid_FOUND
             COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
             COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_boxdfm_2d_trias-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_trias_ug-00064.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_trias_ug params.input -Problem.Name test_2p_boxdfm_trias_ug -Grid.File grids/durlofsky.msh")

# tetrahedra alu
dumux_add_test(NAME test_2p_boxdfm_tets_alu
              LABELS porousmediumflow 2p boxdfm
              SOURCES main.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<2,3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_boxdfm_3d_tets-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_tets_alu-00013.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_tets_alu params.input -Problem.Name test_2p_boxdfm_tets_alu -Grid.File grids/tets_3d.msh -TimeLoop.TEnd 450")

# tetrahedra ug
dumux_add_test(NAME test_2p_boxdfm_tets_ug
              LABELS porousmediumflow 2p boxdfm
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<3>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<2,3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_boxdfm_3d_tets-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_tets_ug-00013.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_boxdfm_tets_ug params.input -Problem.Name test_2p_boxdfm_tets_ug -Grid.File grids/tets_3d.msh -TimeLoop.TEnd 450")
