# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

# Case A draft: planar phase-field evaporation with Stefan-limit QoI logging.
dumux_add_test(NAME test_ff_stokes_2p_stefan
               SOURCES main.cc
               LABELS freeflow navierstokes
               TIMEOUT 3600
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_2p_stefan
               CMD_ARGS params.input)
