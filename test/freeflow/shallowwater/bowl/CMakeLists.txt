# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input plot_results.py)

dumux_add_test(NAME test_shallowwater_bowl
               SOURCES main.cc
               LABELS shallowwater
               COMPILE_DEFINITIONS TYPETAG=Bowl LINEARSOLVER=AMGBiCGSTABIstlSolver
               TIMEOUT 3600
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_bowl-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/bowl-00013.vtu
                          --zeroThreshold {"velocityY":1e-14}
                          --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_bowl")

dumux_add_test(NAME test_shallowwater_bowl_parallel
               TARGET test_shallowwater_bowl
               LABELS shallowwater
               TIMEOUT 3600
               CMAKE_GUARD MPI_FOUND
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_bowl-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/s0002-bowl-parallel-00013.pvtu
                          --zeroThreshold {"velocityY":1e-14} --ignore "process rank"
                          --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_bowl -Problem.Name bowl-parallel")

# the second-order flux damps the oscillation less: the tolerance for the L1 error of the water depth with
# respect to the exact solution after two periods lies between the errors of the first- and second-order flux
dumux_add_test(NAME test_shallowwater_bowl_secondorder
               SOURCES main.cc
               LABELS shallowwater
               COMPILE_DEFINITIONS TYPETAG=BowlSecondOrder LINEARSOLVER=ILUBiCGSTABIstlSolver
               CMD_ARGS params.input -Problem.Name bowl_secondorder -Problem.ErrorTolerance 1e-3)

# BDF2 damps the oscillation less than implicit Euler: the tolerance lies between the errors of both schemes
# with the second-order flux
dumux_add_test(NAME test_shallowwater_bowl_secondorder_bdf2
               TARGET test_shallowwater_bowl_secondorder
               LABELS shallowwater
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_bowl_secondorder
               CMD_ARGS params.input -Problem.Name bowl_secondorder_bdf2 -TimeLoop.Scheme BDF2 -Problem.ErrorTolerance 7e-4)
