dune_symlink_to_source_files(FILES params.input convergencetest.py)
dumux_add_test(NAME test_ff_stokes_channel_pipe
               LABELS freeflow navierstokes
               SOURCES main.cc
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ./convergencetest.py
               CMD_ARGS test_ff_stokes_channel_pipe)
