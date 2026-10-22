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
              MPI_RANKS 2
              TIMEOUT 1200
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_2p_adaptive_tpfa_parallel-00001.pvtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_tpfa params.input -Problem.Name test_2p_adaptive_tpfa_parallel")

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
              MPI_RANKS 2
              TIMEOUT 1200
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_pointsource_adaptive_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_2p_pointsource_adaptive_tpfa_parallel-00001.pvtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_pointsource_adaptive_tpfa params.input -Problem.Name test_2p_pointsource_adaptive_tpfa_parallel")

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
# The box vertex output of the fields that depend on the material law (capillary pressure,
# nonwetting-phase pressure, mobilities) is not unique on the lens boundary: such a vertex has
# sub-control volumes in both materials and the value of the last visited element is written.
# The element order depends on the grid implementation and on the domain decomposition, so these
# fields are not compared.
dumux_add_test(NAME test_2p_adaptive_box_simplex
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveBoxSimplex
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --ignore p_napl pc mob_aq mob_napl
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_box_simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box_simplex-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box_simplex params.input -Grid.CellType Simplex -Problem.Name test_2p_adaptive_box_simplex")

# box with conforming simplex refinement, parallel
dumux_add_test(NAME test_2p_adaptive_box_simplex_parallel
              TARGET test_2p_adaptive_box_simplex
              LABELS porousmediumflow 2p parallel
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleAdaptiveBoxSimplex
              CMAKE_GUARD "( dune-alugrid_FOUND AND MPI_FOUND )"
              MPI_RANKS 2
              TIMEOUT 1200
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy --zeroThreshold {"process rank":100}
                       --ignore p_napl pc mob_aq mob_napl
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_adaptive_box_simplex-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_2p_adaptive_box_simplex_parallel-00001.pvtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_adaptive_box_simplex params.input -Grid.CellType Simplex -Problem.Name test_2p_adaptive_box_simplex_parallel")
