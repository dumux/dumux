add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py")

# executable for tpfa tests
add_executable(test_1p_convergence_analytic_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_1p_convergence_analytic_tpfa PUBLIC "TYPETAG=OnePConvergenceTpfa" PUBLIC "GRIDTYPE=Dune::UGGrid<2>")

# using tpfa and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_tpfa_structured
              TARGET test_1p_convergence_analytic_tpfa
              CMAKE_GUARD "( dune-functions_FOUND )"
              LABELS porousmediumflow 1p
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_convergence_analytic_tpfa params.input
                       -Problem.Name test_1p_convergence_analytic_tpfa_structured)
