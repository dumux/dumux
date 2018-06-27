add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dune_add_test(NAME test_gridcreator_cake
              SOURCES test_gridcreator_cake.cc
              CMAKE_GUARD "( dune-uggrid_FOUND OR dune-alugrid_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_cake"
                       --files ${CMAKE_SOURCE_DIR}/test/references/cake-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/cake-00000.vtu)

dune_add_test(NAME test_gridcreator_subgrid
              SOURCES test_gridcreator_subgrid.cc
              CMAKE_GUARD dune-subgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_subgrid"
                       --files ${CMAKE_SOURCE_DIR}/test/references/subgrid-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/subgrid_three.vtu)

dune_symlink_to_source_files(FILES "test_gridcreator_gmsh_3d.input"
                                   "test_gridcreator_gmsh_e_markers.input"
                                   "test_gridcreator_dgf_e_markers.input")

add_executable(test_gridcreator_gmsh_3d_alu EXCLUDE_FROM_ALL test_gridcreator_gmsh_3d.cc)
target_compile_definitions(test_gridcreator_gmsh_3d_alu PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming>)

dune_add_test(NAME test_gridcreator_gmsh_3d_alu_sequential
              TARGET test_gridcreator_gmsh_3d_alu
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_3d_alu -Problem.Name bifurcation_alu -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_alu-00000.vtu)

dune_add_test(NAME test_gridcreator_gmsh_3d_alu_parallel
              TARGET test_gridcreator_gmsh_3d_alu
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_3d_alu -Problem.Name bifurcation_alu_parallel -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-bifurcation_alu_parallel-00000.pvtu)

add_executable(test_gridcreator_gmsh_3d_ug EXCLUDE_FROM_ALL test_gridcreator_gmsh_3d.cc)
target_compile_definitions(test_gridcreator_gmsh_3d_ug PUBLIC GRIDTYPE=Dune::UGGrid<3>)

dune_add_test(NAME test_gridcreator_gmsh_3d_ug_sequential
              TARGET test_gridcreator_gmsh_3d_ug
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_3d_ug -Problem.Name bifurcation_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation_ug-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_ug-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/bifurcation_ug-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_ug-00001.vtu)

dune_add_test(NAME test_gridcreator_gmsh_3d_ug_parallel
              TARGET test_gridcreator_gmsh_3d_ug
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_3d_ug -Problem.Name bifurcation_ug_parallel"
                       --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation_ug-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-bifurcation_ug_parallel-00000.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/bifurcation_ug-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-bifurcation_ug_parallel-00001.pvtu)

add_executable(test_gridcreator_gmsh_e_markers_alu EXCLUDE_FROM_ALL test_gridcreator_gmsh_e_markers.cc)
target_compile_definitions(test_gridcreator_gmsh_e_markers_alu PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>)

dune_add_test(NAME test_gridcreator_gmsh_e_markers_alu_sequential
              TARGET test_gridcreator_gmsh_e_markers_alu
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_alu -Problem.Name fracture_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00001.vtu)

dune_add_test(NAME test_gridcreator_gmsh_e_markers_alu_parallel
              TARGET test_gridcreator_gmsh_e_markers_alu
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_alu -Problem.Name fracture_alu_parallel"
                      --files ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_alu_parallel-00000.pvtu
                              ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference-refined.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_alu_parallel-00001.pvtu)

add_executable(test_gridcreator_gmsh_e_markers_ug EXCLUDE_FROM_ALL test_gridcreator_gmsh_e_markers.cc)
target_compile_definitions(test_gridcreator_gmsh_e_markers_ug PUBLIC GRIDTYPE=Dune::UGGrid<2>)

dune_add_test(NAME test_gridcreator_gmsh_e_markers_ug_sequential
              TARGET test_gridcreator_gmsh_e_markers_ug
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_ug -Problem.Name fracture_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00001.vtu)

dune_add_test(NAME test_gridcreator_gmsh_e_markers_ug_parallel
              TARGET test_gridcreator_gmsh_e_markers_ug
              CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_ug -Problem.Name fracture_ug_parallel"
                      --files ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_ug_parallel-00000.pvtu
                              ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference-refined.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_ug_parallel-00001.pvtu)

add_executable(test_gridcreator_dgf_e_markers_alu EXCLUDE_FROM_ALL test_gridcreator_dgf_e_markers.cc)
target_compile_definitions(test_gridcreator_dgf_e_markers_alu PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>)

dune_add_test(NAME test_gridcreator_dgf_e_markers_alu_sequential
              TARGET test_gridcreator_dgf_e_markers_alu
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_dgf_e_markers_alu -Problem.Name co2_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-00001.vtu)

dune_add_test(NAME test_gridcreator_dgf_e_markers_alu_parallel
              TARGET test_gridcreator_dgf_e_markers_alu
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_dgf_e_markers_alu -Problem.Name co2_alu_parallel"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-00000.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-00001.pvtu)

add_executable(test_gridcreator_dgf_e_markers_ug EXCLUDE_FROM_ALL test_gridcreator_dgf_e_markers.cc)
target_compile_definitions(test_gridcreator_dgf_e_markers_ug PUBLIC GRIDTYPE=Dune::UGGrid<2>)

dune_add_test(NAME test_gridcreator_dgf_e_markers_ug_sequential
              TARGET test_gridcreator_dgf_e_markers_ug
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_dgf_e_markers_ug -Problem.Name co2_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_ug-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_ug-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-00001.vtu)

dune_add_test(NAME test_gridcreator_dgf_e_markers_ug_parallel
              TARGET test_gridcreator_dgf_e_markers_ug
              CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"rank":100}
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_dgf_e_markers_ug -Problem.Name co2_ug_parallel"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_ug-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-00000.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_ug-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-00001.pvtu)
