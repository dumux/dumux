add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

dumux_add_test(NAME test_1p_convergence_structured
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
              COMPILE_DEFINITIONS TYPETAG=OnePConvergence
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_convergence_structured params.input
                       -Problem.Name test_1p_convergence_structured)

dumux_add_test(NAME test_1p_convergence_unstructured
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2>
              COMPILE_DEFINITIONS TYPETAG=OnePConvergence
              TIMEOUT 1000
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_convergence_unstructured params.input
                       -Problem.Name test_1p_convergence_unstructured
                       -Grid.File ../../incompressible/grids/randomlydistorted.dgf)
