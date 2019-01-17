dune_add_test(NAME test_ff_navierstokes_sincos_steady
              SOURCES main.cc
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_steady-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_steady-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_steady params.input
                             -Problem.Name test_ff_navierstokes_sincos_steady")

dune_symlink_to_source_files(FILES "params.input")
