add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

# executable for tpfa tests
add_executable(test_1p_convergence_analytic_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_convergence_analytic_tpfa PUBLIC "TYPETAG=OnePConvergenceTpfa" PUBLIC "GRIDTYPE=Dune::UGGrid<2>")

# executable for mpfa tests
add_executable(test_1p_convergence_analytic_mpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_convergence_analytic_mpfa PUBLIC "TYPETAG=OnePConvergenceMpfa" PUBLIC "GRIDTYPE=Dune::UGGrid<2>")

# using tpfa and structured grid with zero off-diagonal permeability entries
dumux_add_test(NAME test_1p_convergence_analytic_tpfa_structured
              TARGET test_1p_convergence_analytic_tpfa
              LABELS porousmediumflow 1p
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_convergence_analytic_tpfa params.input
                       -Problem.Name test_1p_convergence_analytic_tpfa_structured
                       -Problem.C 0.0)

# using mpfa and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_mpfa_structured
              TARGET test_1p_convergence_analytic_mpfa
              LABELS porousmediumflow 1p
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_convergence_analytic_mpfa params.input
                       -Problem.Name test_1p_convergence_analytic_mpfa_structured)

# using mpfa and unstructured grid
dumux_add_test(NAME test_1p_convergence_analytic_mpfa_unstructured
              TARGET test_1p_convergence_analytic_mpfa
              LABELS porousmediumflow 1p
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_convergence_analytic_mpfa params.input
                       -Problem.Name test_1p_convergence_analytic_mpfa_unstructured
                       -Grid.File ../../incompressible/grids/randomlydistorted.dgf)
