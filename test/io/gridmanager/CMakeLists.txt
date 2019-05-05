add_input_file_links()
dune_symlink_to_source_files(FILES grids)

add_executable(test_gridmanager_cake EXCLUDE_FROM_ALL test_gridmanager_cake.cc)

dune_add_test(NAME test_gridmanager_cake_360
              TARGET test_gridmanager_cake
              LABELS unit io
              CMAKE_GUARD "( dune-uggrid_FOUND OR dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake -Grid.Name 360"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_360-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/cake-3d-360.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_360-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/cake-2d-360.vtu)

dune_add_test(NAME test_gridmanager_cake_210
              TARGET test_gridmanager_cake
              LABELS unit io
              CMAKE_GUARD "( dune-uggrid_FOUND OR dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake -Grid.Name 210 -Grid.Angular1 '0.0 210.0'"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_210-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/cake-3d-210.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_210-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/cake-2d-210.vtu)

dune_add_test(NAME test_gridmanager_subgrid
              SOURCES test_gridmanager_subgrid.cc
              LABELS unit io
              CMAKE_GUARD dune-subgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_subgrid"
                       --files ${CMAKE_SOURCE_DIR}/test/references/subgrid-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/subgrid_three.vtu)

dune_symlink_to_source_files(FILES "test_gridmanager_gmsh_3d.input"
                                   "test_gridmanager_gmsh_e_markers.input"
                                   "test_gridmanager_dgf.input")

add_executable(test_gridmanager_gmsh_3d_alu EXCLUDE_FROM_ALL test_gridmanager_gmsh_3d.cc)
target_compile_definitions(test_gridmanager_gmsh_3d_alu PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming>)

dune_add_test(NAME test_gridmanager_gmsh_3d_alu_sequential
              TARGET test_gridmanager_gmsh_3d_alu
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_alu -Problem.Name bifurcation_alu -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_alu-00000.vtu)

dune_add_test(NAME test_gridmanager_gmsh_3d_alu_parallel
              TARGET test_gridmanager_gmsh_3d_alu
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_alu -Problem.Name bifurcation_alu_parallel -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-bifurcation_alu_parallel-00000.pvtu)

add_executable(test_gridmanager_gmsh_3d_ug EXCLUDE_FROM_ALL test_gridmanager_gmsh_3d.cc)
target_compile_definitions(test_gridmanager_gmsh_3d_ug PUBLIC GRIDTYPE=Dune::UGGrid<3>)

dune_add_test(NAME test_gridmanager_gmsh_3d_ug_sequential
              TARGET test_gridmanager_gmsh_3d_ug
              LABELS unit io
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_ug -Problem.Name bifurcation_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_ug-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-bifurcation-3d-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_ug-00001.vtu)

dune_add_test(NAME test_gridmanager_gmsh_3d_ug_parallel
              TARGET test_gridmanager_gmsh_3d_ug
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_ug -Problem.Name bifurcation_ug_parallel -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-bifurcation_ug_parallel-00000.pvtu)

add_executable(test_gridmanager_gmsh_e_markers_alu EXCLUDE_FROM_ALL test_gridmanager_gmsh_e_markers.cc)
target_compile_definitions(test_gridmanager_gmsh_e_markers_alu PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>)

dune_add_test(NAME test_gridmanager_gmsh_e_markers_alu_sequential
              TARGET test_gridmanager_gmsh_e_markers_alu
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_alu -Problem.Name fracture_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00001.vtu)

dune_add_test(NAME test_gridmanager_gmsh_e_markers_alu_parallel
              TARGET test_gridmanager_gmsh_e_markers_alu
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_alu -Problem.Name fracture_alu_parallel"
                      --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_alu_parallel-00000.pvtu
                              ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference-refined.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_alu_parallel-00001.pvtu)

add_executable(test_gridmanager_gmsh_e_markers_ug EXCLUDE_FROM_ALL test_gridmanager_gmsh_e_markers.cc)
target_compile_definitions(test_gridmanager_gmsh_e_markers_ug PUBLIC GRIDTYPE=Dune::UGGrid<2>)

dune_add_test(NAME test_gridmanager_gmsh_e_markers_ug_sequential
              TARGET test_gridmanager_gmsh_e_markers_ug
              LABELS unit io
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_ug -Problem.Name fracture_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00001.vtu)

dune_add_test(NAME test_gridmanager_gmsh_e_markers_ug_parallel
              TARGET test_gridmanager_gmsh_e_markers_ug
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_ug -Problem.Name fracture_ug_parallel"
                      --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_ug_parallel-00000.pvtu
                              ${CMAKE_SOURCE_DIR}/test/references/gridmanager-fracture-reference-refined.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_ug_parallel-00001.pvtu)

add_executable(test_gridmanager_dgf_alu EXCLUDE_FROM_ALL test_gridmanager_dgf.cc)
target_compile_definitions(test_gridmanager_dgf_alu PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>)

dune_add_test(NAME test_gridmanager_dgf_alu_sequential
              TARGET test_gridmanager_dgf_alu
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_alu -Problem.Name co2_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-simpl-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-element-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-simpl-element-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-element-00001.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-simpl-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-vertex-00000.vtu)

dune_add_test(NAME test_gridmanager_dgf_alu_parallel
              TARGET test_gridmanager_dgf_alu
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_alu -Problem.Name co2_alu_parallel"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-simpl-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-element-00000.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-simpl-element-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-element-00001.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-simpl-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-vertex-00000.pvtu)

add_executable(test_gridmanager_dgf_ug EXCLUDE_FROM_ALL test_gridmanager_dgf.cc)
target_compile_definitions(test_gridmanager_dgf_ug PUBLIC GRIDTYPE=Dune::UGGrid<2>)

dune_add_test(NAME test_gridmanager_dgf_ug_sequential
              TARGET test_gridmanager_dgf_ug
              LABELS unit io
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_ug -Problem.Name co2_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-quad-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-element-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-quad-element-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-element-00001.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-quad-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-vertex-00000.vtu)

dune_add_test(NAME test_gridmanager_dgf_ug_parallel
              TARGET test_gridmanager_dgf_ug
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_ug -Problem.Name co2_ug_parallel"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-quad-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-element-00000.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-quad-element-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-element-00001.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridmanager-co2-quad-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-vertex-00000.pvtu)

dune_add_test(NAME test_gridmanager_vtk_foam
              SOURCES test_gridmanager_vtk.cc
              LABELS unit io
              COMPILE_DEFINITIONS FOAMGRID
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_vtk_foam
                       --files ${CMAKE_CURRENT_BINARY_DIR}/grids/lupine.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test-gridmanager-vtk-foamgrid-0.vtp)

dune_add_test(NAME test_gridmanager_vtk_ug
              SOURCES test_gridmanager_vtk.cc
              LABELS unit io
              COMPILE_DEFINITIONS UGGRID
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_vtk_ug
                       --files ${CMAKE_CURRENT_BINARY_DIR}/grids/lens.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test-gridmanager-vtk-uggrid-0.vtu)
