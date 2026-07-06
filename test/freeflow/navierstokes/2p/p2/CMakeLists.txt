# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input rect.msh)

# P2-P1-P2-P2 Taylor-Hood Cahn-Hilliard / Navier-Stokes rising bubble.
# Build+run target (no regression reference yet; the scheme is under development).
dumux_add_test(NAME test_ff_stokes_2p_p2_taylorhood
               SOURCES main.cc
               LABELS freeflow navierstokes
               TIMEOUT 3600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_2p_p2_taylorhood
               CMD_ARGS params.input)
