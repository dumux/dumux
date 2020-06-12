dune_symlink_to_source_files(FILES params.input)
dumux_add_test(NAME test_ff_stokes_channel_pipe
               LABELS freeflow navierstokes
               SOURCES main.cc
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pipe-00000.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pipe-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_pipe params.input"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-12})
