dune_symlink_to_source_files(FILES params.input cylinder_channel.msh cylinder_channel_quad.msh)

# Navier-Stokes version of the test (Re=20)
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary
               SOURCES main.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND AND HAVE_MPI )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_stationary.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary -Problem.EnableInertiaTerms true")

# Stokes version of the test (no inertia)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary
               TARGET test_ff_navierstokes_dfg_benchmark_stationary
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND AND HAVE_MPI )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary -Problem.EnableInertiaTerms false -LinearSolver.UseIterativeSolver true")
