# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES network.dgf params.input params_gridmanager.input params_neumann.input)

dumux_add_test(NAME test_1p_network1d3d_parallel_tpfa
               LABELS porousmediumflow 1p parallel
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkCCTpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params.input)

dumux_add_test(NAME test_1p_network1d3d_parallel_tpfa_bisection
               LABELS porousmediumflow 1p parallel
               TARGET test_1p_network1d3d_parallel_tpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 2 4
               TIMEOUT 600
               COMMAND ./test_1p_network1d3d_parallel_tpfa
               CMD_ARGS params.input -Grid.PartitionByBulk false -Grid.Overlap 2)

dumux_add_test(NAME test_1p_network1d3d_parallel_box
               LABELS porousmediumflow 1p parallel
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkBox
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params.input)

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

dumux_add_test(NAME test_1p_network1d3d_gridmanager
               LABELS porousmediumflow 1p parallel
               SOURCES main_gridmanager.cc
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               CMD_ARGS params_gridmanager.input)

dumux_add_test(NAME test_1p_network1d3d_parallel_tpfa_refined
               LABELS porousmediumflow 1p parallel
               TARGET test_1p_network1d3d_parallel_tpfa
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               COMMAND ./test_1p_network1d3d_parallel_tpfa
               CMD_ARGS params.input -Grid.MaxElementLength 2e-6)

dumux_add_test(NAME test_1p_network1d3d_parallel_box_refined
               LABELS porousmediumflow 1p parallel
               TARGET test_1p_network1d3d_parallel_box
               CMAKE_GUARD "( dune-foamgrid_FOUND AND MPI_FOUND )"
               MPI_RANKS 1 2 4
               TIMEOUT 600
               COMMAND ./test_1p_network1d3d_parallel_box
               CMD_ARGS params.input -Grid.MaxElementLength 2e-6)
