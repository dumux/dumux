add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

add_executable(test_md_boundary_darcy1p_freeflow1p_convtest EXCLUDE_FROM_ALL main.cc)

dune_add_test(NAME test_md_boundary_darcy1p_stokes1p_convtest
              TARGET test_md_boundary_darcy1p_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1p_freeflow1p_convtest params.input
                        -Problem.TestCase ShiueExampleTwo
                        -Darcy.SpatialParams.Permeability 1.0)

dune_add_test(NAME test_md_boundary_darcy1p_navierstokes1p_convtest
              TARGET test_md_boundary_darcy1p_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1p_freeflow1p_convtest params.input
                        -Problem.TestCase Schneider
                        -FreeFlow.Problem.EnableInertiaTerms true
                        -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true)
