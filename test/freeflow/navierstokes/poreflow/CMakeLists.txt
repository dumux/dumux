add_input_file_links()
dune_symlink_to_source_files(FILES poreflow.pbm)

dumux_add_test(NAME test_ff_poreflow
               SOURCES main.cc
               LABELS freeflow navierstokes subgrid
               CMAKE_GUARD HAVE_UMFPACK
               CMAKE_GUARD dune-subgrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_poreflow-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_poreflow-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_poreflow params.input
                             -Problem.OutletCondition Outflow
                             -Problem.Name test_ff_poreflow")
