# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(momentum)

dune_symlink_to_source_files(FILES params.input convergencetest.py)
dumux_add_test(NAME test_ff_stokes_channel_pipe
               LABELS freeflow navierstokes
               SOURCES main.cc
               LABELS freeflow
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ./convergencetest.py
               CMD_ARGS test_ff_stokes_channel_pipe params.input
                                  -Problem.PrintConvergenceTestFile true)
