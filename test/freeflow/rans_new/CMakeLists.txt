dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_ff_rans_twoeq
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_rans_twoeq-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_twoeq-00009.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_twoeq params.input
                             -Problem.Name test_ff_rans_twoeq")
