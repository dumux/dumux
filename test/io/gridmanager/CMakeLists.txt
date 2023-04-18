# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

add_executable(test_gridmanager_cake_ug EXCLUDE_FROM_ALL test_gridmanager_cake.cc)
target_compile_definitions(test_gridmanager_cake_ug PUBLIC "USEUG=1")

add_executable(test_gridmanager_cake_alu EXCLUDE_FROM_ALL test_gridmanager_cake.cc)

dumux_add_test(NAME test_gridmanager_cake_360_ug
               TARGET test_gridmanager_cake_ug
               LABELS unit io
               CMAKE_GUARD dune-uggrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake_ug -Grid.Name ug-360"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-3d-ug-360.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-2d-ug-360.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_gridmanager_cake_360_ug_parallel
               TARGET test_gridmanager_cake_ug
               LABELS unit io
               CMAKE_GUARD "( MPI_FOUND AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore  "process rank"
                        --command "${MPIEXEC} -np 3 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake_ug -Grid.Name ug-360"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0003-cake-3d-ug-360.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0003-cake-2d-ug-360.pvtu
                        --relative 1e-6)

dumux_add_test(NAME test_gridmanager_cake_360_alu
               TARGET test_gridmanager_cake_alu
               LABELS unit io
               CMAKE_GUARD dune-alugrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake_alu -Grid.Name alu-360"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-3d-alu-360.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-2d-alu-360.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_gridmanager_cake_360_alu_parallel
               TARGET test_gridmanager_cake_alu
               LABELS unit io
               CMAKE_GUARD "( MPI_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy --ignore "process rank"
                        --command "${MPIEXEC} -np 3 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake_alu -Grid.Name alu-360"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0003-cake-3d-alu-360.pvtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_360-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0003-cake-2d-alu-360.pvtu
                        --relative 1e-6)

dumux_add_test(NAME test_gridmanager_cake_210_ug
               TARGET test_gridmanager_cake_ug
               LABELS unit io
               CMAKE_GUARD dune-uggrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake_ug -Grid.Name ug-210 -Grid.Angular1 '0.0 210.0'"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_210-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-3d-ug-210.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_210-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-2d-ug-210.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_gridmanager_cake_210_alu
               TARGET test_gridmanager_cake_alu
               LABELS unit io
               CMAKE_GUARD dune-alugrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake_alu -Grid.Name alu-210 -Grid.Angular1 '0.0 210.0'"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_210-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-3d-alu-210.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_210-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-2d-alu-210.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_gridmanager_cake_360_nowell
               TARGET test_gridmanager_cake_ug
               LABELS unit io
               CMAKE_GUARD dune-uggrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_cake_ug -Grid.Name nowell-360 -Grid.Radial0 '0.0 0.3 1.0'"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_3d_360_nowell-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-3d-nowell-360.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_cake_2d_360_nowell-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/cake-2d-nowell-360.vtu
                        --relative 1e-6)

dumux_add_test(NAME test_gridmanager_subgrid
              SOURCES test_gridmanager_subgrid.cc
              LABELS unit io
              CMAKE_GUARD dune-subgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_subgrid"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_subgrid-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/subgrid_circle_yasp.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager_subgrid_binary_image-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/subgrid_binary_image.vtu)

dune_symlink_to_source_files(FILES "test_gridmanager_gmsh_3d.input"
                                   "test_gridmanager_gmsh_e_markers.input"
                                   "test_gridmanager_dgf.input")

add_executable(test_gridmanager_gmsh_3d_alu EXCLUDE_FROM_ALL test_gridmanager_gmsh_3d.cc)
target_compile_definitions(test_gridmanager_gmsh_3d_alu PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming>)

dumux_add_test(NAME test_gridmanager_gmsh_3d_alu_sequential
              TARGET test_gridmanager_gmsh_3d_alu
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_alu -Problem.Name bifurcation_alu -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_alu-00000.vtu)

dumux_add_test(NAME test_gridmanager_gmsh_3d_alu_parallel
              TARGET test_gridmanager_gmsh_3d_alu
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --ignore "rank"
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_alu -Problem.Name bifurcation_alu_parallel -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-bifurcation_alu_parallel-00000.pvtu)

add_executable(test_gridmanager_gmsh_3d_ug EXCLUDE_FROM_ALL test_gridmanager_gmsh_3d.cc)
target_compile_definitions(test_gridmanager_gmsh_3d_ug PUBLIC GRIDTYPE=Dune::UGGrid<3>)

dumux_add_test(NAME test_gridmanager_gmsh_3d_ug_sequential
              TARGET test_gridmanager_gmsh_3d_ug
              LABELS unit io
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_ug -Problem.Name bifurcation_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_ug-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-bifurcation-3d-refined-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/bifurcation_ug-00001.vtu)

dumux_add_test(NAME test_gridmanager_gmsh_3d_ug_parallel
              TARGET test_gridmanager_gmsh_3d_ug
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --ignore "rank"
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_3d_ug -Problem.Name bifurcation_ug_parallel -Grid.Refine false"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-bifurcation-3d-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-bifurcation_ug_parallel-00000.pvtu)

add_executable(test_gridmanager_gmsh_e_markers_alu EXCLUDE_FROM_ALL test_gridmanager_gmsh_e_markers.cc)
target_compile_definitions(test_gridmanager_gmsh_e_markers_alu PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>)

dumux_add_test(NAME test_gridmanager_gmsh_e_markers_alu_sequential
              TARGET test_gridmanager_gmsh_e_markers_alu
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_alu -Problem.Name fracture_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-refined-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_alu-00001.vtu)

dumux_add_test(NAME test_gridmanager_gmsh_e_markers_alu_parallel
              TARGET test_gridmanager_gmsh_e_markers_alu
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --ignore "rank"
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_alu -Problem.Name fracture_alu_parallel"
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_alu_parallel-00000.pvtu
                              ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-refined-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_alu_parallel-00001.pvtu)

add_executable(test_gridmanager_gmsh_e_markers_ug EXCLUDE_FROM_ALL test_gridmanager_gmsh_e_markers.cc)
target_compile_definitions(test_gridmanager_gmsh_e_markers_ug PUBLIC GRIDTYPE=Dune::UGGrid<2>)

dumux_add_test(NAME test_gridmanager_gmsh_e_markers_ug_sequential
              TARGET test_gridmanager_gmsh_e_markers_ug
              LABELS unit io
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_ug -Problem.Name fracture_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-refined-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/fracture_ug-00001.vtu)

dumux_add_test(NAME test_gridmanager_gmsh_e_markers_ug_parallel
              TARGET test_gridmanager_gmsh_e_markers_ug
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --ignore "rank"
                      --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_gmsh_e_markers_ug -Problem.Name fracture_ug_parallel"
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_ug_parallel-00000.pvtu
                              ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-fracture-refined-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/s0002-fracture_ug_parallel-00001.pvtu)

add_executable(test_gridmanager_dgf_alu EXCLUDE_FROM_ALL test_gridmanager_dgf.cc)
target_compile_definitions(test_gridmanager_dgf_alu PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>)

dumux_add_test(NAME test_gridmanager_dgf_alu_sequential
              TARGET test_gridmanager_dgf_alu
              LABELS unit io
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_alu -Problem.Name co2_alu"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-simpl-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-element-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-simpl-element-refined-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-element-00001.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-simpl-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_alu-vertex-00000.vtu)

dumux_add_test(NAME test_gridmanager_dgf_alu_parallel
              TARGET test_gridmanager_dgf_alu
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --ignore "rank"
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_alu -Problem.Name co2_alu_parallel"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-simpl-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-element-00000.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-simpl-element-refined-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-element-00001.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-simpl-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_alu_parallel-vertex-00000.pvtu)

add_executable(test_gridmanager_dgf_ug EXCLUDE_FROM_ALL test_gridmanager_dgf.cc)
target_compile_definitions(test_gridmanager_dgf_ug PUBLIC GRIDTYPE=Dune::UGGrid<2>)

dumux_add_test(NAME test_gridmanager_dgf_ug_sequential
              TARGET test_gridmanager_dgf_ug
              LABELS unit io
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_ug -Problem.Name co2_ug"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-quad-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-element-00000.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-quad-element-refined-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-element-00001.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-quad-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/co2_ug-vertex-00000.vtu)

dumux_add_test(NAME test_gridmanager_dgf_ug_parallel
              TARGET test_gridmanager_dgf_ug
              LABELS unit io parallel
              TIMEOUT 300
              CMAKE_GUARD "( dune-uggrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --ignore "rank"
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_dgf_ug -Problem.Name co2_ug_parallel"
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-quad-element-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-element-00000.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-quad-element-refined-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-element-00001.pvtu
                               ${CMAKE_SOURCE_DIR}/test/references/test_gridmanager-co2-quad-vertex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-co2_ug_parallel-vertex-00000.pvtu)

dumux_add_test(NAME test_gridmanager_vtk_foam
              SOURCES test_gridmanager_vtk.cc
              LABELS unit io
              COMPILE_DEFINITIONS FOAMGRID
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_vtk_foam
                       --files ${CMAKE_CURRENT_BINARY_DIR}/grids/lupine.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test-gridmanager-vtk-foamgrid-0.vtp)

dumux_add_test(NAME test_gridmanager_vtk_ug
              SOURCES test_gridmanager_vtk.cc
              LABELS unit io
              COMPILE_DEFINITIONS UGGRID
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --command ${CMAKE_CURRENT_BINARY_DIR}/test_gridmanager_vtk_ug
                       --files ${CMAKE_CURRENT_BINARY_DIR}/grids/lens.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test-gridmanager-vtk-uggrid-0.vtu)


dumux_add_test(NAME test_gmshboundaryflag
              SOURCES test_gmshboundaryflag.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false
              LABELS unit io
              CMAKE_GUARD "( dune-alugrid_FOUND )"
              )

dumux_add_test(NAME test_gmshboundaryflag_caching
              SOURCES test_gmshboundaryflag.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true
              LABELS unit io
              CMAKE_GUARD "( dune-alugrid_FOUND )"
              )

dumux_add_test(NAME test_gridmanager_mmesh
              SOURCES test_gridmanager_mmesh.cc
              LABELS unit io
              CMAKE_GUARD "( dune-mmesh_FOUND )"
              )

add_executable(test_pnm_gridmanager EXCLUDE_FROM_ALL test_pnm_gridmanager.cc)

dumux_add_test(NAME test_pnm_gridmanager_constant_params
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-3d-constant-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/generic-3dgrid-constant-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-3d-constant-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/dgf-3dgrid-constant-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-2d-constant-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/generic-2dgrid-constant-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-2d-constant-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/dgf-2dgrid-constant-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-1d-constant-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/generic-1dgrid-constant-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-1d-constant-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/dgf-1dgrid-constant-params.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager test_pnm_gridmanager-constant-params.input")

dumux_add_test(NAME test_pnm_gridmanager_random_params
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-3d-random-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/generic-3dgrid-random-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-3d-random-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/dgf-3dgrid-random-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-2d-random-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/generic-2dgrid-random-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-2d-random-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/dgf-2dgrid-random-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-1d-random-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/generic-1dgrid-random-params.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_pnm_gridmanager-1d-random-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/dgf-1dgrid-random-params.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager test_pnm_gridmanager-random-params.input")

dumux_add_test(NAME test_pnm_gridmanager_remove_throats_on_boundary_none
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager
               CMD_ARGS test_remove-troats-on-boundary.input)

dumux_add_test(NAME test_pnm_gridmanager_remove_throats_on_boundary_0
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager
               CMD_ARGS test_remove-troats-on-boundary.input -Grid.RemoveThroatsOnBoundary 0)

dumux_add_test(NAME test_pnm_gridmanager_remove_throats_on_boundary_1
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager
               CMD_ARGS test_remove-troats-on-boundary.input -Grid.RemoveThroatsOnBoundary 1)

dumux_add_test(NAME test_pnm_gridmanager_remove_throats_on_boundary_2
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager
               CMD_ARGS test_remove-troats-on-boundary.input -Grid.RemoveThroatsOnBoundary 2)

dumux_add_test(NAME test_pnm_gridmanager_remove_throats_on_boundary_3
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager
               CMD_ARGS test_remove-troats-on-boundary.input -Grid.RemoveThroatsOnBoundary 3)

dumux_add_test(NAME test_pnm_gridmanager_remove_throats_on_boundary_all
               LABELS unit io porenetwork
               TARGET test_pnm_gridmanager
               CMAKE_GUARD dune-foamgrid_FOUND
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_gridmanager
               CMD_ARGS test_remove-troats-on-boundary.input -Grid.RemoveThroatsOnBoundary "0 1 2 3")

dumux_add_test(SOURCES test_pnm_subgriddata.cc
               LABELS unit io porenetwork
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-subgrid_FOUND )")
