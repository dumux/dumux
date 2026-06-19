# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES network.dgf params.input params_gridmanager.input params_neumann.input)

# 1d network (Dune::FoamGrid<1,3>) distributed over a coarse YaspGrid: the test self-checks that
# the parallel (overlapping) solution equals the sequential one, so it only needs the exit code.
dumux_add_test(NAME test_1p_network1d3d_parallel_tpfa
               LABELS porousmediumflow 1p parallel
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkCCTpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params.input)

dumux_add_test(NAME test_1p_network1d3d_parallel_box
               LABELS porousmediumflow 1p parallel
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkBox
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params.input)

# mixed Neumann/Dirichlet BC on a distributed connected path: parallel == sequential
dumux_add_test(NAME test_1p_network1d3d_neumann_tpfa
               LABELS porousmediumflow 1p parallel
               SOURCES main_neumann.cc
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkCCTpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params_neumann.input)
dumux_add_test(NAME test_1p_network1d3d_neumann_box
               LABELS porousmediumflow 1p parallel
               SOURCES main_neumann.cc
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkBox
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params_neumann.input)

# drives the distribution through the standard Dumux::GridManager (Grid.ParallelPartitioning)
dumux_add_test(NAME test_1p_network1d3d_gridmanager
               LABELS porousmediumflow 1p parallel
               SOURCES main_gridmanager.cc
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params_gridmanager.input)
