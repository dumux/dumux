dune_add_test(NAME test_ff_navierstokes_kovasznay
              SOURCES main.cc
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_kovasznay-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_kovasznay params.input
                             -Problem.Name test_ff_navierstokes_kovasznay")
dune_symlink_to_source_files(FILES "params.input")
