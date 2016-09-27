add_input_file_links()
dune_symlink_to_source_files(FILES grids)

add_executable(test_gridcreator_gmsh EXCLUDE_FROM_ALL test_gridcreator_gmsh.cc)

dune_add_test(NAME test_gridcreator_gmsh
              TARGET test_gridcreator_gmsh
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh"
                       --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation-00001.vtu)

dune_add_test(NAME test_gridcreator_gmsh_parallel
              TARGET test_gridcreator_gmsh
              CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference-parallel.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-bifurcation-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/bifurcation-reference-refined-parallel.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0000-bifurcation-00001.vtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh")

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

dune_symlink_to_source_files(FILES "test_gridcreator_gmsh_e_markers.input"
                                   "test_gridcreator_dgf_e_markers.input")

# UG currently doesn't implement element data load balancing
# add_executable(test_gridcreator_gmsh_e_markers_ug EXCLUDE_FROM_ALL test_gridcreator_gmsh_e_markers_ug.cc)

# dune_add_test(NAME test_gridcreator_gmsh_e_markers_ug
              # TARGET test_gridcreator_gmsh_e_markers_ug
              # CMAKE_GUARD dune-uggrid_FOUND
              # COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              # CMD_ARGS --script fuzzy
                      # --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_ug"
                      # --files ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference.vtu
                              # ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00000.vtu
                              # ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference-refined.vtu
                            # ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00001.vtu)
# 
# dune_add_test(NAME test_gridcreator_gmsh_e_markers_ug_parallel
              # TARGET test_gridcreator_gmsh_e_markers_ug
              # CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              # COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              # CMD_ARGS --script fuzzy
                       # --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_ug"
                       # --files ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference-parallel.vtu
                               # ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-fracture_ug-00000.vtu
                               # ${CMAKE_SOURCE_DIR}/test/references/fracture_ug-reference-refined-parallel.vtu
                             # ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-fracture_ug-00001.vtu)

add_executable(test_gridcreator_gmsh_e_markers_alu EXCLUDE_FROM_ALL test_gridcreator_gmsh_e_markers_alu.cc)

dune_add_test(NAME test_gridcreator_gmsh_e_markers_alu
              TARGET test_gridcreator_gmsh_e_markers_alu
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00001.vtu)

dune_add_test(NAME test_gridcreator_gmsh_e_markers_alu_parallel
              TARGET test_gridcreator_gmsh_e_markers_alu
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_gmsh_e_markers_alu"
                      --files ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference-parallel.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-fracture_alu-00000.vtu
                              ${CMAKE_SOURCE_DIR}/test/references/fracture_alu-reference-refined-parallel.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-fracture_alu-00001.vtu)

add_executable(test_gridcreator_dgf_e_markers_alu EXCLUDE_FROM_ALL test_gridcreator_dgf_e_markers_alu.cc)

dune_add_test(NAME test_gridcreator_dgf_e_markers_alu
              TARGET test_gridcreator_dgf_e_markers_alu
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_dgf_e_markers_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-00001.vtu)

dune_add_test(NAME test_gridcreator_dgf_e_markers_alu_parallel
              TARGET test_gridcreator_dgf_e_markers_alu
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridcreator_dgf_e_markers_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference-parallel.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-co2_alu-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/gridcreator-co2_alu-reference-parallel-refined.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-p0001-co2_alu-00001.vtu)
