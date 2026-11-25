# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_md_embedded_1d3d_loadbalance_alugrid
               LABELS multidomain multidomain_embedded parallel
               SOURCES main.cc
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 300
               CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_loadbalance_alugrid_refined
               TARGET test_md_embedded_1d3d_loadbalance_alugrid
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 300
               CMD_ARGS params.input -Grid.Refinement 1)

dumux_add_test(NAME test_md_embedded_1d3d_loadbalance_scotch
               TARGET test_md_embedded_1d3d_loadbalance_alugrid
               LABELS multidomain multidomain_embedded parallel
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND MPI_FOUND AND PTSCOTCH_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 300
               CMD_ARGS params.input -LoadBalance.Method Scotch)
