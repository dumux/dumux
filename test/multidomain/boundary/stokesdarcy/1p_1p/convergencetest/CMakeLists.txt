add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

add_executable(test_md_boundary_darcy1ptpfa_freeflow1p_convtest EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy1ptpfa_freeflow1p_convtest PUBLIC DARCYTYPETAG=DarcyOnePCC)

add_executable(test_md_boundary_darcy1pbox_freeflow1p_convtest EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_boundary_darcy1pbox_freeflow1p_convtest PUBLIC DARCYTYPETAG=DarcyOnePBox)

dune_add_test(NAME test_md_boundary_darcy1ptpfa_stokes1p_convtest
              TARGET test_md_boundary_darcy1ptpfa_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1ptpfa_freeflow1p_convtest params.input
                        -Problem.TestCase ShiueExampleTwo
                        -Darcy.SpatialParams.Permeability 1.0)

dune_add_test(NAME test_md_boundary_darcy1ptpfa_navierstokes1p_convtest
              TARGET test_md_boundary_darcy1ptpfa_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1ptpfa_freeflow1p_convtest params.input
                        -Problem.TestCase Schneider
                        -FreeFlow.Problem.EnableInertiaTerms true
                        -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true)

dune_add_test(NAME test_md_boundary_darcy1pbox_stokes1p_convtest
              TARGET test_md_boundary_darcy1pbox_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1pbox_freeflow1p_convtest params.input
                          -Problem.TestCase ShiueExampleTwo
                          -Darcy.SpatialParams.Permeability 1.0)

dune_add_test(NAME test_md_boundary_darcy1pbox_navierstokes1p_convtest
              TARGET test_md_boundary_darcy1pbox_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1pbox_freeflow1p_convtest params.input
                          -Problem.TestCase Schneider
                          -FreeFlow.Problem.EnableInertiaTerms true
                          -FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph true)

dune_add_test(NAME test_md_boundary_darcy1pbox_stokes1p_convtest_bj
              TARGET test_md_boundary_darcy1pbox_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1pbox_freeflow1p_convtest params.input
                          -Problem.TestCase BJSymmetrized
                          -Problem.SlipCondition BJ
                          -Darcy.SpatialParams.Porosity 1.0
                          -Darcy.SpatialParams.Permeability "1.0 1.0 2.0")

dune_add_test(NAME test_md_boundary_darcy1p_stokes1p_newic
              TARGET test_md_boundary_darcy1pbox_freeflow1p_convtest
              LABELS multidomain multidomain_boundary stokesdarcy
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS  test_md_boundary_darcy1pbox_freeflow1p_convtest params.input
                          -Problem.TestCase NewICNonSymmetrized
                          -Problem.SlipCondition ER
                          -Problem.NewIc true
                          -FreeFlow.EnableUnsymmetrizedVelocityGradientForIC true
                          -FreeFlow.EnableUnsymmetrizedVelocityGradient true
                          -Darcy.InterfaceParams.EpsInterface 1.0
                          -Darcy.InterfaceParams.N_s_bl 0.5
                          -Darcy.InterfaceParams.N_1_bl -3.0
                          -Darcy.InterfaceParams.M_bl "1.0 1.0"
                          -Darcy.SpatialParams.Porosity 1.0
                          -Darcy.SpatialParams.Permeability "1.0 1.0 0.0"
                          -Component.LiquidKinematicViscosity 2.0)