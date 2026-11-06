# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

add_executable(test_md_boundary_ff1p_pm1p_convtest EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq2fe_box EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ1Bubble
                                  TYPETAG_FFMASS=FreeFlowOnePMassBox
                                  TYPETAG_PMMASS=DarcyOnePBoxNewInterface)

target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ1BubbleHybrid
                                  TYPETAG_FFMASS=FreeFlowOnePMassBox
                                  TYPETAG_PMMASS=DarcyOnePBoxNewInterface)

target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ2Hybrid
                                  TYPETAG_FFMASS=FreeFlowOnePMassBox
                                  TYPETAG_PMMASS=DarcyOnePBoxNewInterface
                                  QUADRATURE_RULE_FFMASS=Dumux::QuadratureRules::DuneQuadrature<4>)

# pure finite element velocity: PQ2 with a Box pressure, and the full Taylor-Hood pair
target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq2fe_box
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ2FE
                                  TYPETAG_FFMASS=FreeFlowOnePMassBox
                                  TYPETAG_PMMASS=DarcyOnePBoxNewInterface
                                  QUADRATURE_RULE_FFMASS=Dumux::QuadratureRules::DuneQuadrature<4>)

target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ2FE
                                  TYPETAG_FFMASS=FreeFlowOnePMassPQ1FE
                                  TYPETAG_PMMASS=DarcyOnePPQ1FE)

# PQ3 velocity and PQ2 pressure, with a PQ2 or a PQ3 Darcy pressure
target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ3Hybrid
                                  TYPETAG_FFMASS=FreeFlowOnePMassPQ2Hybrid
                                  TYPETAG_PMMASS=DarcyOnePPQ2Hybrid)

target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ3Hybrid
                                  TYPETAG_FFMASS=FreeFlowOnePMassPQ2Hybrid
                                  TYPETAG_PMMASS=DarcyOnePPQ3Hybrid)

target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  TYPETAG_FFMOMENTUM=FreeFlowOnePMomentumPQ3FE
                                  TYPETAG_FFMASS=FreeFlowOnePMassPQ2FE
                                  TYPETAG_PMMASS=DarcyOnePPQ3FE)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.Name freeflow_staggered -Darcy.Problem.Name darcy_staggered"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_shiue1-errors.csv
                       --check "L2dof(p_ff) all 1.7 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.Name freeflow_pq1bubble_box -Darcy.Problem.Name darcy_pq1bubble_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue1-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.Name freeflow_pq1bubble_hybrid_box -Darcy.Problem.Name darcy_pq1bubble_hybrid_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue1-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '20 20' -FreeFlow.Grid.Cells '20 20' -FreeFlow.Problem.Name freeflow_pq2_hybrid_box -Darcy.Problem.Name darcy_pq2_hybrid_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue1-errors.csv
                       --check "H1(v) all 1.9 2.1"
                       --check "L2(p_ff) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2fe_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2fe_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq2fe_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq2fe_box_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '20 20' -FreeFlow.Grid.Cells '20 20' -FreeFlow.Problem.Name freeflow_pq2fe_box -Darcy.Problem.Name darcy_pq2fe_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq2fe_box_shiue1-errors.csv
                       --check "H1(v) all 1.9 2.1"
                       --check "L2(p_ff) all 1.9 2.2"
                       --check "L2(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '20 20' -FreeFlow.Grid.Cells '20 20' -FreeFlow.Problem.Name freeflow_taylorhood_fe -Darcy.Problem.Name darcy_taylorhood_fe"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe_shiue1-errors.csv
                       --check "H1(v) all 1.9 2.1"
                       --check "L2(p_ff) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 3000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '10 10' -FreeFlow.Grid.Cells '10 10' -FreeFlow.Problem.Name freeflow_pq3_pq2_hybrid -Darcy.Problem.Name darcy_pq3_pq2_hybrid"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid_shiue1-errors.csv
                       --check "H1(v) all 2.7 3.1"
                       --check "L2(p_ff) all 2.7 3.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 3000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '10 10' -FreeFlow.Grid.Cells '10 10' -FreeFlow.Problem.Name freeflow_pq3_pq3_hybrid -Darcy.Problem.Name darcy_pq3_pq3_hybrid"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid_shiue1-errors.csv
                       --check "H1(v) all 2.9 3.1"
                       --check "L2(p_ff) all 2.9 3.1"
                       --check "L2(p_pm) all 3.8 4.2")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 3000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe_shiue1 -Problem.TestCase ShiueExampleOne -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '10 10' -FreeFlow.Grid.Cells '10 10' -FreeFlow.Problem.Name freeflow_pq3_pq3_fe -Darcy.Problem.Name darcy_pq3_pq3_fe"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe_shiue1-errors.csv
                       --check "H1(v) all 2.9 3.1"
                       --check "L2(p_ff) all 2.9 3.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_shiue2 -Problem.TestCase ShiueExampleTwo -Darcy.Problem.BoundaryConditions Mixed -FreeFlow.Problem.BoundaryConditions Mixed -FreeFlow.Problem.Name freeflow_staggered -Darcy.Problem.Name darcy_staggered"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_shiue2-errors.csv
                       --check "L2dof(p_ff) all 1.7 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue2 -Problem.TestCase ShiueExampleTwo -Darcy.Problem.BoundaryConditions Mixed -FreeFlow.Problem.BoundaryConditions Mixed -FreeFlow.Problem.Name freeflow_pq1bubble_box -Darcy.Problem.Name darcy_pq1bubble_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue2-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue2 -Problem.TestCase ShiueExampleTwo -Darcy.Problem.BoundaryConditions Mixed -FreeFlow.Problem.BoundaryConditions Mixed -FreeFlow.Problem.Name freeflow_pq1bubble_hybrid_box -Darcy.Problem.Name darcy_pq1bubble_hybrid_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue2-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue2 -Problem.TestCase ShiueExampleTwo -Darcy.Problem.BoundaryConditions Mixed -FreeFlow.Problem.BoundaryConditions Mixed -Darcy.Grid.Cells '20 20' -FreeFlow.Grid.Cells '20 20' -FreeFlow.Problem.Name freeflow_pq2_hybrid_box -Darcy.Problem.Name darcy_pq2_hybrid_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue2-errors.csv
                       --check "L2(p_ff) all 1.9 2.1"
                       --check "L2dof(p_ff) all 1.9 2.1"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 1.9 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_rybak
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_rybak -Problem.TestCase Rybak -Darcy.Problem.BoundaryConditions Neumann -FreeFlow.Problem.BoundaryConditions Stress -FreeFlow.Problem.Name freeflow_staggered -Darcy.Problem.Name darcy_staggered"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_rybak-errors.csv
                       --check "L2dof(p_ff) all 1.9 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak -Problem.TestCase Rybak -Darcy.Problem.BoundaryConditions Neumann -FreeFlow.Problem.BoundaryConditions Stress -FreeFlow.Problem.Name freeflow_pq1bubble_box -Darcy.Problem.Name darcy_pq1bubble_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak -Problem.TestCase Rybak -Darcy.Problem.BoundaryConditions Neumann -FreeFlow.Problem.BoundaryConditions Stress -Darcy.Grid.Cells '20 20' -FreeFlow.Grid.Cells '20 20' -FreeFlow.Problem.Name freeflow_pq2_hybrid_box -Darcy.Problem.Name darcy_pq2_hybrid_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak-errors.csv
                       --check "L2(p_ff) all 1.9 2.1"
                       --check "L2dof(p_ff) all 1.9 2.1"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 1.9 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_navierstokes
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_navierstokes -Problem.TestCase Schneider -FreeFlow.Problem.EnableInertiaTerms true -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true -Darcy.Problem.BoundaryConditions Mixed -FreeFlow.Problem.BoundaryConditions Mixed -FreeFlow.Problem.Name freeflow_staggered -Darcy.Problem.Name darcy_staggered"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_navierstokes-errors.csv
                       --check "L2dof(p_ff) all 1.9 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_navierstokes
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_navierstokes -Problem.TestCase Schneider -FreeFlow.Problem.EnableInertiaTerms true -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true -Darcy.Problem.BoundaryConditions Mixed -FreeFlow.Problem.BoundaryConditions Mixed -FreeFlow.Problem.Name freeflow_pq1bubble_box -Darcy.Problem.Name darcy_pq1bubble_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_navierstokes-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_navierstokes
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_navierstokes -Problem.TestCase Schneider -FreeFlow.Problem.EnableInertiaTerms true -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true -Darcy.Problem.BoundaryConditions Mixed -FreeFlow.Problem.BoundaryConditions Mixed -Darcy.Grid.Cells '20 20' -FreeFlow.Grid.Cells '20 20' -FreeFlow.Problem.Name freeflow_pq2_hybrid_box -Darcy.Problem.Name darcy_pq2_hybrid_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_navierstokes-errors.csv
                       --check "L2(p_ff) all 1.9 2.1"
                       --check "L2dof(p_ff) all 1.9 2.1"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 1.9 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_beaversjoseph
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_beaversjoseph -Problem.TestCase BeaversJoseph -Darcy.SpatialParams.Permeability 1e-2 -FreeFlow.Problem.SlipCondition BJ -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.Name freeflow_pq1bubble_box -Darcy.Problem.Name darcy_pq1bubble_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_beaversjoseph-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

# the exact solution is contained in the discrete spaces
dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_beaversjoseph
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_beaversjoseph -Problem.TestCase BeaversJoseph -Darcy.SpatialParams.Permeability 1e-2 -FreeFlow.Problem.SlipCondition BJ -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '10 10' -FreeFlow.Grid.Cells '10 10' -FreeFlow.Problem.Name freeflow_pq2_hybrid_box -Darcy.Problem.Name darcy_pq2_hybrid_box"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_beaversjoseph-errors.csv
                       --round-off 1e-12
                       --check "L2(p_ff) all 1.9 2.1"
                       --check "L2dof(p_ff) all 1.9 2.1"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 1.9 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak_bj
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak_bj -Problem.TestCase Rybak -FreeFlow.Problem.SlipCondition BJ -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.Name freeflow_pq1bubble_box_bj -Darcy.Problem.Name darcy_pq1bubble_box_bj"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak_bj-errors.csv
                       --check "L2(p_ff) all 1.4 1.7"
                       --check "L2dof(p_ff) all 1.4 1.7"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 0.9 1.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak_bj
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input -Problem.Name test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak_bj -Problem.TestCase Rybak -FreeFlow.Problem.SlipCondition BJ -Darcy.Problem.BoundaryConditions Dirichlet -FreeFlow.Problem.BoundaryConditions Dirichlet -Darcy.Grid.Cells '20 20' -FreeFlow.Grid.Cells '20 20' -FreeFlow.Problem.Name freeflow_pq2_hybrid_box_bj -Darcy.Problem.Name darcy_pq2_hybrid_box_bj"
                       --levels 0 1 2
                       --level-args "-Grid.Refinement {}"
                       --errors test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak_bj-errors.csv
                       --check "L2(p_ff) all 1.9 2.1"
                       --check "L2dof(p_ff) all 1.9 2.1"
                       --check "L2(v) all 1.9 2.1"
                       --check "H1(v) all 1.9 2.1"
                       --check "L2dof(v) all 1.9 2.1"
                       --check "L2(p_pm) all 1.9 2.1"
                       --check "L2dof(p_pm) all 1.9 2.1")
