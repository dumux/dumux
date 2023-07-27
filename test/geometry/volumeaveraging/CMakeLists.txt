add_input_file_links()
dune_symlink_to_source_files(FILES test.pbm)


dumux_add_test(NAME test_volumeaveraging
               SOURCES main.cc
               LABELS unit volumeaveraging
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_volumeaveraging_base-ref.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/BaseGrid-00000.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/test_volumeaveraging_extended-ref.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/ExtendedGrid-00000.vtu
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_volumeaveraging params.input")

dumux_add_test(NAME test_volumeaveraging_subgrid
               SOURCES main_subgrid.cc
               LABELS unit volumeaveraging
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/Averaging_SourceGrid-ref.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/SourceGrid-00000.vtu
                                 ${CMAKE_SOURCE_DIR}/test/references/Averaging_TargetGrid-ref.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/TargetGrid-00000.vtu
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_volumeaveraging_subgrid params_subgrid.input
                         -Source.Grid.Cells '150 150'
                         -Target.Grid.Cells '50 50' ")
