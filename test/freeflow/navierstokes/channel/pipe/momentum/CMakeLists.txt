dune_symlink_to_source_files(FILES params.input convergencetest.py)
dumux_add_test(NAME test_ff_stokes_channel_pipe_momenum_diamond
               LABELS freeflow navierstokes
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=PipeFlowDiamond
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ./convergencetest.py
               CMD_ARGS test_ff_stokes_channel_pipe_momenum_diamond)

dumux_add_test(NAME test_ff_stokes_channel_pipe_momenum_cvfe
               LABELS freeflow navierstokes
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=PipeFlowCvfe
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ./convergencetest.py
               CMD_ARGS test_ff_stokes_channel_pipe_momenum_cvfe)