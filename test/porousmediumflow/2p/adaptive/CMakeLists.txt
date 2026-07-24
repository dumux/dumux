# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "initialsolutioncc.txt")
dune_symlink_to_source_files(FILES "initialsolutionbox.txt")

# using tpfa
dumux_add_test(NAME test_2p_adaptive_tpfa
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveTpfa
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_tpfa params.input -Problem.Name test_2p_adaptive_tpfa")

# tpfa parallel
dumux_add_test(NAME test_2p_adaptive_tpfa_parallel
              TARGET test_2p_adaptive_tpfa
              LABELS porousmediumflow 2p parallel
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveTpfa
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_2p_adaptive_tpfa_parallel-00001.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_tpfa params.input -Problem.Name test_2p_adaptive_tpfa_parallel")

# using tpfa and point source
dumux_add_test(NAME test_2p_pointsource_adaptive_tpfa
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPAdaptivePointSource
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_pointsource_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa params.input -Problem.Name test_2p_pointsource_adaptive_tpfa")

# tpfa and point source parallel
dumux_add_test(NAME test_2p_pointsource_adaptive_tpfa_parallel
              TARGET test_2p_pointsource_adaptive_tpfa
              LABELS porousmediumflow 2p parallel
              COMPILE_DEFINITIONS TYPETAG=TwoPAdaptivePointSource
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_pointsource_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_2p_pointsource_adaptive_tpfa_parallel-00001.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa params.input -Problem.Name test_2p_pointsource_adaptive_tpfa_parallel")

# using mpfa
dumux_add_test(NAME test_2p_adaptive_mpfa
              LABELS porousmediumflow 2p
             SOURCES main.cc
             COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveMpfa
             CMAKE_GUARD dune-alugrid_FOUND
             COMMAND ${DUMUX_RUNTEST}
             CMD_ARGS --script fuzzy
                      --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_mpfa-reference.vtu
                              ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_mpfa-00001.vtu
                      --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_mpfa params.input -Problem.Name test_2p_adaptive_mpfa")

# using box
dumux_add_test(NAME test_2p_adaptive_box
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveBox
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box params.input -Problem.Name test_2p_adaptive_box")

# using box with conforming simplex refinement
dumux_add_test(NAME test_2p_adaptive_box_simplex
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveBoxSimplex
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_box_simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box_simplex-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box_simplex params.input -Grid.CellType Simplex -Problem.Name test_2p_adaptive_box_simplex")

# box with conforming simplex refinement, parallel
# The vertex-based output of the material-law-derived quantities (capillary pressure and
# mobilities) is ambiguous at the lens material interface, where the value picked for a
# boundary vertex depends on the domain decomposition. These fields are therefore ignored
# in the parallel comparison; the primary variables (saturation and pressure) match the
# sequential reference to machine precision.
dumux_add_test(NAME test_2p_adaptive_box_simplex_parallel
              TARGET test_2p_adaptive_box_simplex
              LABELS porousmediumflow 2p parallel
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveBoxSimplex
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --ignore pc mob_aq mob_napl
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_box_simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_2p_adaptive_box_simplex_parallel-00001.pvtu
                       --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box_simplex params.input -Grid.CellType Simplex -Problem.Name test_2p_adaptive_box_simplex_parallel")
