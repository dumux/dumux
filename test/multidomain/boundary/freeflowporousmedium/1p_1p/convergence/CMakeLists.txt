# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

add_executable(test_md_boundary_ff1p_pm1p_convtest EXCLUDE_FROM_ALL main.cc)

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
