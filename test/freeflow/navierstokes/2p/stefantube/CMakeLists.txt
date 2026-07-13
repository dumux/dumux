# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_3d.input params_3d_filament.input)

# Stefan-tube slit-capillary evaporation with curved meniscus and reduced QoIs.
dumux_add_test(NAME test_ff_stokes_2p_stefantube
               SOURCES main.cc
               LABELS freeflow navierstokes
               TIMEOUT 3600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_2p_stefantube
               CMD_ARGS params.input)

dumux_add_test(NAME test_ff_stokes_2p_stefantube_3d
               SOURCES main.cc
               COMPILE_DEFINITIONS STEFANTUBE_DIM=3
               LABELS freeflow navierstokes
               TIMEOUT 3600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_2p_stefantube_3d
               CMD_ARGS params_3d.input)
