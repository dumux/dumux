dune_symlink_to_source_files(FILES "params.input" "channel.msh" "channel_quad.msh")
dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform
               SOURCES main.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform params.input
                             -Problem.Name test_ff_stokes_channel_3d_nonuniform")

dumux_add_test(NAME test_ff_stokes_channel_3d_nonuniform_simplex
               TARGET test_ff_stokes_channel_3d_nonuniform
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_3d_nonuniform_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform_simplex-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_3d_nonuniform params.input
                             -Problem.Name test_ff_stokes_channel_3d_nonuniform_simplex
                             -Grid.File channel.msh")
