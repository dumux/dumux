dune_add_test(NAME test_ff_navierstokes_sincos_unsteady
              SOURCES main.cc
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_unsteady-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_unsteady-00017.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_unsteady params.input
                             -Problem.Name test_ff_navierstokes_sincos_unsteady")

dune_symlink_to_source_files(FILES "params.input")
