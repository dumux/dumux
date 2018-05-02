dune_symlink_to_source_files(FILES "test_2p.input" "grids")

# quadrilaterals yasp
dune_add_test(NAME test_2pboxdfm_quads_yasp
              SOURCES test_2pboxdfm.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/2pboxdfm_2d_quads-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/2pboxdfm_quads_yasp-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pboxdfm_quads_yasp test_2p.input -Problem.Name 2pboxdfm_quads_yasp")

# quadrilaterals alu grid
dune_add_test(NAME test_2pboxdfm_quads_alu
              SOURCES test_2pboxdfm.cc
              CMAKE_GUARD dune-alugrid_FOUND
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/2pboxdfm_2d_quads_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/2pboxdfm_quads_alu-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pboxdfm_quads_alu test_2p.input -Problem.Name 2pboxdfm_quads_alu")

# quadrilaterals ug grid
dune_add_test(NAME test_2pboxdfm_quads_ug
              SOURCES test_2pboxdfm.cc
              CMAKE_GUARD dune-uggrid_FOUND
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/2pboxdfm_2d_quads-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/2pboxdfm_quads_ug-00011.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pboxdfm_quads_ug test_2p.input -Problem.Name 2pboxdfm_quads_ug")

# triangles alu
dune_add_test(NAME test_2pboxdfm_trias_alu
             SOURCES test_2pboxdfm.cc
             CMAKE_GUARD dune-foamgrid_FOUND
             CMAKE_GUARD dune-alugrid_FOUND
             COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
             COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/2pboxdfm_2d_trias_alu-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/2pboxdfm_trias_alu-00007.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pboxdfm_trias_alu test_2p.input -Problem.Name 2pboxdfm_trias_alu -Grid.File grids/2d_trias.msh")

# triangles ug
dune_add_test(NAME test_2pboxdfm_trias_ug
             SOURCES test_2pboxdfm.cc
             CMAKE_GUARD dune-uggrid_FOUND
             CMAKE_GUARD dune-foamgrid_FOUND
             COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
             COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<1,2>
             COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/2pboxdfm_2d_trias_ug-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/2pboxdfm_trias_ug-00007.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pboxdfm_trias_ug test_2p.input -Problem.Name 2pboxdfm_trias_ug -Grid.File grids/2d_trias.msh")

# tetrahedra alu
dune_add_test(NAME test_2pboxdfm_tets_alu
              SOURCES test_2pboxdfm.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              CMAKE_GUARD dune-alugrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<2,3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/2pboxdfm_2d_tets_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/2pboxdfm_tets_alu-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pboxdfm_tets_alu test_2p.input -Problem.Name 2pboxdfm_tets_alu -Grid.File grids/3d_simplex.msh")

# tetrahedra ug
dune_add_test(NAME test_2pboxdfm_tets_ug
              SOURCES test_2pboxdfm.cc
              CMAKE_GUARD dune-uggrid_FOUND
              CMAKE_GUARD dune-foamgrid_FOUND
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<3>
              COMPILE_DEFINITIONS FRACTUREGRIDTYPE=Dune::FoamGrid<2,3>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/2pboxdfm_2d_tets_ug-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/2pboxdfm_tets_ug-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pboxdfm_tets_ug test_2p.input -Problem.Name 2pboxdfm_tets_ug -Grid.File grids/3d_simplex.msh")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
test_2p_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/incompressible)
