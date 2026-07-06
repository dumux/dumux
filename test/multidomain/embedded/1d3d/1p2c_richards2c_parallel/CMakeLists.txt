# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# MPI-parallel transient root-soil benchmark (1p2c Richards / 1pnc root) on a realistic DGF
# network, overlapping (Strategy A) decomposition solved with the parallel multidomain MUMPS
# backend. Compares the distributed transient solution (and migrated radii) against a replicated
# sequential reference and returns 0 only if they match (see main.cc).

dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c_parallel
              LABELS multidomain multidomain_embedded richards 1pnc parallel mumps
              SOURCES main.cc
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p2c_richards2c_parallel
              CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c_parallel_2p
              TARGET test_md_embedded_1d3d_1p2c_richards2c_parallel
              LABELS multidomain multidomain_embedded richards 1pnc parallel mumps
              CMAKE_GUARD "( MPI_FOUND AND MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ${MPIEXEC_EXECUTABLE}
              CMD_ARGS ${MPIEXEC_NUMPROC_FLAG} 2
                       ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p2c_richards2c_parallel params.input)

# Same benchmark, but the distributed solve uses a block-diagonal (Restricted Additive Schwarz)
# preconditioner: each rank factorizes only its local block with MUMPS inside a parallel BiCGSTAB.
# The built-in comparison against the replicated sequential reference validates correctness.
dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c_parallel_ras
              LABELS multidomain multidomain_embedded richards 1pnc parallel mumps
              SOURCES main.cc
              COMPILE_DEFINITIONS USE_RAS_SOLVER
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p2c_richards2c_parallel_ras
              CMD_ARGS params.input)

dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c_parallel_ras_2p
              TARGET test_md_embedded_1d3d_1p2c_richards2c_parallel_ras
              LABELS multidomain multidomain_embedded richards 1pnc parallel mumps
              CMAKE_GUARD "( MPI_FOUND AND MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ${MPIEXEC_EXECUTABLE}
              CMD_ARGS ${MPIEXEC_NUMPROC_FLAG} 2
                       ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p2c_richards2c_parallel_ras params.input)

# production driver: distributed run over the real time span with VTK output, mass balance and
# timing (used for the scaling study; also a short smoke test here)
dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c_production
              LABELS multidomain multidomain_embedded richards 1pnc parallel mumps
              SOURCES main_production.cc
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p2c_richards2c_production
              CMD_ARGS params_production.input -TimeLoop.TEnd 17280 -Problem.EnableVtkOutput false
                       -Soil.Grid.Cells "10 10 10")

# production driver with the block-diagonal (RAS) solver instead of the global MUMPS direct solve
dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c_production_ras
              LABELS multidomain multidomain_embedded richards 1pnc parallel mumps
              SOURCES main_production.cc
              COMPILE_DEFINITIONS USE_RAS_SOLVER
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p2c_richards2c_production_ras
              CMD_ARGS params_production.input -TimeLoop.TEnd 17280 -Problem.EnableVtkOutput false
                       -Soil.Grid.Cells "10 10 10")

# production driver with the original sequential field-block-diagonal ILU0 BiCGSTAB solver
# (sequential baseline for comparing the new MUMPS-based solvers against)
dumux_add_test(NAME test_md_embedded_1d3d_1p2c_richards2c_production_blockdiagilu
              LABELS multidomain multidomain_embedded richards 1pnc mumps
              SOURCES main_production.cc
              COMPILE_DEFINITIONS USE_BLOCKDIAGILU_SOLVER
              CMAKE_GUARD "( MUMPS_FOUND AND dune-foamgrid_FOUND )"
              COMMAND ./test_md_embedded_1d3d_1p2c_richards2c_production_blockdiagilu
              CMD_ARGS params_production.input -TimeLoop.TEnd 17280 -Problem.EnableVtkOutput false
                       -Soil.Grid.Cells "10 10 10")

dune_symlink_to_source_files(FILES "params.input" "params_production.input")
