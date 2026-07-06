# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# MPI-parallel 1d-3d embedded coupling (overlapping decomposition) solved with the parallel
# multidomain MUMPS direct solver. The test compares the distributed solution against a
# replicated sequential reference and returns 0 only if they match (see main.cc).

# Milestone 1: line coupling (cc-tpfa, no circle average / extended stencil).
# Sequential smoke run (1 rank): exercises the same code path; parallel == reference trivially.
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_parallel_line
              LABELS multidomain multidomain_embedded 1p parallel mumps
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p1p_parallel_line
              CMD_ARGS params.input)

# The real parallel check on 2 ranks.
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_parallel_line_2p
              TARGET test_md_embedded_1d3d_1p1p_parallel_line
              LABELS multidomain multidomain_embedded 1p parallel mumps
              CMAKE_GUARD "( MPI_FOUND AND MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ${MPIEXEC_EXECUTABLE}
              CMD_ARGS ${MPIEXEC_NUMPROC_FLAG} 2
                       ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_parallel_line params.input)

# Milestone 2: average coupling (circle average + extended bulk-bulk source stencil).
# The bulk is split across the circle plane (x) with an overlap wide enough to cover the
# vessel radius (radius 0.3, h = 0.1 -> overlap >= 3), so the circle average of an owned
# boundary cell reaches its neighbours locally.
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_parallel_average
              LABELS multidomain multidomain_embedded 1p parallel mumps
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p1p_parallel_average
              CMD_ARGS params.input -Tissue.Grid.Overlap 3 -MixedDimension.CheckCircleOverlap true)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_parallel_average_2p
              TARGET test_md_embedded_1d3d_1p1p_parallel_average
              LABELS multidomain multidomain_embedded 1p parallel mumps
              CMAKE_GUARD "( MPI_FOUND AND MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ${MPIEXEC_EXECUTABLE}
              CMD_ARGS ${MPIEXEC_NUMPROC_FLAG} 2
                       ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_parallel_average params.input
                       -Tissue.Grid.Overlap 3 -Tissue.Grid.Partitioning "2 1 1"
                       -MixedDimension.CheckCircleOverlap true)

# Realistic DGF network (lupine root system) with per-element radius migrated by global id.
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_parallel_network
              LABELS multidomain multidomain_embedded 1p parallel mumps
              SOURCES main_network.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=NetworkCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p1p_parallel_network
              CMD_ARGS params_network.input)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_parallel_network_2p
              TARGET test_md_embedded_1d3d_1p1p_parallel_network
              LABELS multidomain multidomain_embedded 1p parallel mumps
              CMAKE_GUARD "( MPI_FOUND AND MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ${MPIEXEC_EXECUTABLE}
              CMD_ARGS ${MPIEXEC_NUMPROC_FLAG} 2
                       ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_parallel_network params_network.input)

dune_symlink_to_source_files(FILES "params.input" "params_network.input")
