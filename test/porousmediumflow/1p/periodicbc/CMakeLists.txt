# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "periodic.dgf")

# test using periodic boundary conditions
dumux_add_test(NAME test_1p_periodic_tpfa
              SOURCES main.cc
              LABELS porousmediumflow 1p
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa FVGEOMCACHING=false
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_tpfa_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa params.input -Problem.Name test_1p_periodic_tpfa")

dumux_add_test(NAME test_1p_periodic_tpfa_parallel
              TARGET test_1p_periodic_tpfa
              LABELS porousmediumflow 1p parallel
              CMAKE_GUARD "( MPI_FOUND AND dune-spgrid_FOUND )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --ignore "process rank"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_tpfa_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_1p_periodic_tpfa_parallel-00001.pvtu
                        --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa params.input -Problem.Name test_1p_periodic_tpfa_parallel")

dumux_add_test(NAME test_1p_periodic_tpfa_caching
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa FVGEOMCACHING=true
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_tpfa_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa_caching-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa params.input -Problem.Name test_1p_periodic_tpfa_caching")

dumux_add_test(NAME test_1p_periodic_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleBox
              CMAKE_GUARD dune-spgrid_FOUND
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_box params.input -Problem.Name test_1p_periodic_box")

dumux_add_test(NAME test_1p_periodic_tpfa_subgrid
              SOURCES main.cc
              LABELS porousmediumflow 1p
              COMPILE_DEFINITIONS TYPETAG=OnePIncompressibleTpfa FVGEOMCACHING=false USESUBGRID=true
              CMAKE_GUARD "( dune-spgrid_FOUND AND dune-subgrid_FOUND )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_tpfa_periodic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa_subgrid-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_periodic_tpfa_subgrid params.input -Problem.Name test_1p_periodic_tpfa_subgrid")
