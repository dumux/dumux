dune_symlink_to_source_files(FILES "convergencetest.py" "params.input" "sphere.msh")

# using box and structured grid
dumux_add_test(NAME test_1p_convergence_analytic_sphere_nltpfa
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=OnePConvergenceWMpfa
               COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<3>
               CMAKE_GUARD dune-uggrid_FOUND
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_analytic_sphere_nltpfa)
