# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "convergencetest.py")

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
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest params.input
                        -Problem.TestCase ShiueExampleOne
                        -Rates "freeflow:L2dof(p) 1.7 2.1" "freeflow:L2dof(v) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.Name freeflow_staggered
                        -Darcy.Problem.Name darcy_staggered)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_box
                        -Darcy.Problem.Name darcy_pq1bubble_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_hybrid_box
                        -Darcy.Problem.Name darcy_pq1bubble_hybrid_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "20 20"
                        -FreeFlow.Grid.Cells "20 20"
                        -Rates "H1(v) 1.9 2.1" "freeflow:L2(p) 1.9 2.1" "darcy:L2(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq2_hybrid_box
                        -Darcy.Problem.Name darcy_pq2_hybrid_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2fe_box_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2fe_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq2fe_box params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "20 20"
                        -FreeFlow.Grid.Cells "20 20"
                        -Rates "H1(v) 1.9 2.1" "freeflow:L2(p) 1.9 2.2" "darcy:L2(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq2fe_box
                        -Darcy.Problem.Name darcy_pq2fe_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_taylorhood_fe params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "20 20"
                        -FreeFlow.Grid.Cells "20 20"
                        -Rates "H1(v) 1.9 2.1" "freeflow:L2(p) 1.9 2.1" "darcy:L2(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_taylorhood_fe
                        -Darcy.Problem.Name darcy_taylorhood_fe)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 3000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq3_pq2_hybrid params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "10 10"
                        -FreeFlow.Grid.Cells "10 10"
                        -Rates "H1(v) 2.7 3.1" "freeflow:L2(p) 2.7 3.1"
                        -FreeFlow.Problem.Name freeflow_pq3_pq2_hybrid
                        -Darcy.Problem.Name darcy_pq3_pq2_hybrid)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 3000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_hybrid params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "10 10"
                        -FreeFlow.Grid.Cells "10 10"
                        -Rates "H1(v) 2.9 3.1" "freeflow:L2(p) 2.9 3.1" "darcy:L2(p) 3.8 4.2"
                        -FreeFlow.Problem.Name freeflow_pq3_pq3_hybrid
                        -Darcy.Problem.Name darcy_pq3_pq3_hybrid)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 3000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq3_pq3_fe params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "10 10"
                        -FreeFlow.Grid.Cells "10 10"
                        -Rates "H1(v) 2.9 3.1" "freeflow:L2(p) 2.9 3.1"
                        -FreeFlow.Problem.Name freeflow_pq3_pq3_fe
                        -Darcy.Problem.Name darcy_pq3_pq3_fe)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest params.input
                        -Problem.TestCase ShiueExampleTwo
                        -Rates "freeflow:L2dof(p) 1.7 2.1" "freeflow:L2dof(v) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.Name freeflow_staggered
                        -Darcy.Problem.Name darcy_staggered)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input
                        -Problem.TestCase ShiueExampleTwo
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_box
                        -Darcy.Problem.Name darcy_pq1bubble_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_hybrid_box params.input
                        -Problem.TestCase ShiueExampleTwo
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_hybrid_box
                        -Darcy.Problem.Name darcy_pq1bubble_hybrid_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input
                        -Problem.TestCase ShiueExampleTwo
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed
                        -Darcy.Grid.Cells "20 20"
                        -FreeFlow.Grid.Cells "20 20"
                        -Rates "L2(p) 1.9 2.1" "L2(v) 1.9 2.1" "H1(v) 1.9 2.1" "L2dof(p) 1.9 2.1" "L2dof(v) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq2_hybrid_box
                        -Darcy.Problem.Name darcy_pq2_hybrid_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_rybak
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest params.input
                        -Problem.TestCase Rybak
                        -Darcy.Problem.BoundaryConditions Neumann
                        -FreeFlow.Problem.BoundaryConditions Stress
                        -FreeFlow.Problem.Name freeflow_staggered
                        -Darcy.Problem.Name darcy_staggered)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input
                        -Problem.TestCase Rybak
                        -Darcy.Problem.BoundaryConditions Neumann
                        -FreeFlow.Problem.BoundaryConditions Stress
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_box
                        -Darcy.Problem.Name darcy_pq1bubble_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input
                        -Problem.TestCase Rybak
                        -Darcy.Problem.BoundaryConditions Neumann
                        -FreeFlow.Problem.BoundaryConditions Stress
                        -Darcy.Grid.Cells "20 20"
                        -FreeFlow.Grid.Cells "20 20"
                        -Rates "L2(p) 1.9 2.1" "L2(v) 1.9 2.1" "H1(v) 1.9 2.1" "L2dof(p) 1.9 2.1" "L2dof(v) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq2_hybrid_box
                        -Darcy.Problem.Name darcy_pq2_hybrid_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_navierstokes
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest params.input
                        -Problem.TestCase Schneider
                        -FreeFlow.Problem.EnableInertiaTerms true
                        -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.Name freeflow_staggered
                        -Darcy.Problem.Name darcy_staggered)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_navierstokes
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input
                        -Problem.TestCase Schneider
                        -FreeFlow.Problem.EnableInertiaTerms true
                        -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_box
                        -Darcy.Problem.Name darcy_pq1bubble_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_navierstokes
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input
                        -Problem.TestCase Schneider
                        -FreeFlow.Problem.EnableInertiaTerms true
                        -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed
                        -Darcy.Grid.Cells "20 20"
                        -FreeFlow.Grid.Cells "20 20"
                        -Rates "L2(p) 1.9 2.1" "L2(v) 1.9 2.1" "H1(v) 1.9 2.1" "L2dof(p) 1.9 2.1" "L2dof(v) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq2_hybrid_box
                        -Darcy.Problem.Name darcy_pq2_hybrid_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_beaversjoseph
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input
                        -Problem.TestCase BeaversJoseph
                        -Darcy.SpatialParams.Permeability 1e-2
                        -FreeFlow.Problem.SlipCondition BJ
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_box
                        -Darcy.Problem.Name darcy_pq1bubble_box)

# the exact solution is contained in the discrete spaces
dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_beaversjoseph
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input
                        -Problem.TestCase BeaversJoseph
                        -Darcy.SpatialParams.Permeability 1e-2
                        -FreeFlow.Problem.SlipCondition BJ
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "10 10"
                        -FreeFlow.Grid.Cells "10 10"
                        -Rates "L2(p) 1.9 2.1" "L2(v) 1.9 2.1" "H1(v) 1.9 2.1" "L2dof(p) 1.9 2.1" "L2dof(v) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq2_hybrid_box
                        -Darcy.Problem.Name darcy_pq2_hybrid_box)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box_rybak_bj
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box params.input
                        -Problem.TestCase Rybak
                        -FreeFlow.Problem.SlipCondition BJ
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Rates "freeflow:L2(p) 1.4 1.7" "freeflow:L2dof(p) 1.4 1.7" "L2(v) 1.9 2.1" "H1(v) 0.9 1.1"
                               "freeflow:L2dof(v) 1.9 2.1" "darcy:L2(p) 1.9 2.1" "darcy:L2dof(p) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq1bubble_box_bj
                        -Darcy.Problem.Name darcy_pq1bubble_box_bj)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box_rybak_bj
              TARGET test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest_pq2_hybrid_box params.input
                        -Problem.TestCase Rybak
                        -FreeFlow.Problem.SlipCondition BJ
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet
                        -Darcy.Grid.Cells "20 20"
                        -FreeFlow.Grid.Cells "20 20"
                        -Rates "L2(p) 1.9 2.1" "L2(v) 1.9 2.1" "H1(v) 1.9 2.1" "L2dof(p) 1.9 2.1" "L2dof(v) 1.9 2.1"
                        -FreeFlow.Problem.Name freeflow_pq2_hybrid_box_bj
                        -Darcy.Problem.Name darcy_pq2_hybrid_box_bj)
