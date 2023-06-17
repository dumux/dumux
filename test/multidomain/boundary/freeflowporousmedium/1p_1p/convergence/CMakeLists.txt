# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

add_executable(test_md_boundary_ff1p_pm1p_convtest EXCLUDE_FROM_ALL main.cc)
add_executable(test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_ff1p_pm1p_convtest_pq1bubble_box
                           PUBLIC GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                                  MOMENTUM_DISCRETIZATION_MODEL=PQ1BubbleModel
                                  NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
                                  MASS_DISCRETIZATION_MODEL=BoxModel)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_shiue1
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest params.input
                        -Problem.TestCase ShiueExampleOne
                        -Darcy.Problem.BoundaryConditions Dirichlet
                        -FreeFlow.Problem.BoundaryConditions Dirichlet)

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
                        -Rates 1.4,1.55,1.95,2.05,1.95,2.05,1.95,2.05)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_shiue2
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest params.input
                        -Problem.TestCase ShiueExampleTwo
                        -Darcy.Problem.BoundaryConditions Mixed
                        -FreeFlow.Problem.BoundaryConditions Mixed)

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
                        -Rates 1.4,1.55,1.95,2.05,1.95,2.05,1.95,2.05)

dumux_add_test(NAME test_md_boundary_ff1p_pm1p_convtest_rybak
              TARGET test_md_boundary_ff1p_pm1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_ff1p_pm1p_convtest params.input
                        -Problem.TestCase Rybak
                        -Darcy.Problem.BoundaryConditions Neumann
                        -FreeFlow.Problem.BoundaryConditions Stress)

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
                        -Rates 1.4,1.55,1.95,2.05,1.95,2.05,1.95,2.05)

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
                        -FreeFlow.Problem.BoundaryConditions Mixed)

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
                        -FreeFlow.Problem.BoundaryConditions Stress
                        -Rates 1.4,1.55,1.95,2.05,1.95,2.05,1.95,2.05)