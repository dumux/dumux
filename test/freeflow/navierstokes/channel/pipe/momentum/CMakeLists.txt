# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py)
dumux_add_test(NAME test_ff_stokes_channel_pipe_momentum_diamond
               LABELS freeflow navierstokes
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=PipeFlowDiamond
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ./convergencetest.py
               CMD_ARGS test_ff_stokes_channel_pipe_momentum_diamond)

dumux_add_test(NAME test_ff_stokes_channel_pipe_momentum_pq1bubble
               LABELS freeflow navierstokes
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=PipeFlowPQ1Bubble
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ./convergencetest.py
               CMD_ARGS test_ff_stokes_channel_pipe_momentum_pq1bubble)
