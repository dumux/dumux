dune_symlink_to_source_files(FILES "convergencetest.py")

add_executable(test_ff_navierstokes_sincos EXCLUDE_FROM_ALL main.cc)

dune_add_test(NAME test_ff_navierstokes_sincos_unsteady
              TARGET test_ff_navierstokes_sincos
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_unsteady-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_unsteady-00017.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos params_unsteady.input
                             -Problem.Name test_ff_navierstokes_sincos_unsteady")

dune_add_test(NAME test_ff_navierstokes_sincos_steady
              TARGET test_ff_navierstokes_sincos
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_steady-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_steady-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos params_steady.input
                             -Problem.Name test_ff_navierstokes_sincos_steady")

dune_add_test(NAME test_ff_navierstokes_sincos_steady_convergence
              TARGET test_ff_navierstokes_sincos
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ./convergencetest.py
              CMD_ARGS       test_ff_navierstokes_sincos params_steady.input  -Grid.Cells "150 150"
                             -Problem.Name test_ff_navierstokes_sincos_steady_convergence)

dune_symlink_to_source_files(FILES "params_steady.input")
dune_symlink_to_source_files(FILES "params_unsteady.input")
