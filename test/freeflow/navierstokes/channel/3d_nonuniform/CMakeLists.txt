dune_symlink_to_source_files(FILES "params.input" "channel.msh")
dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform
               SOURCES main.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform params.input
                             -Problem.Name test_ff_stokes_channel_3d_nonuniform"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-12})
