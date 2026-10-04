# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)
dumux_add_test(NAME test_ff_navierstokes_channel_annulus
               LABELS freeflow navierstokes
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_channel_annulus
               CMD_ARGS params.input)

dumux_add_test(NAME test_ff_stokes_channel_annulus
               TARGET test_ff_navierstokes_channel_annulus
               LABELS freeflow navierstokes
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_channel_annulus
               CMD_ARGS params.input -Problem.Name test_ff_stokes_channel_annulus
                        -Problem.EnableInertiaTerms false -Problem.MaxPressureError 2e-3)
