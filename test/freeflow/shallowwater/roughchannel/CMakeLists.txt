# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(limitednikuradse)
dune_symlink_to_source_files(FILES params.input convergencetest.py plot_results.py)

dumux_add_test(NAME test_shallowwater_roughchannel
               SOURCES main.cc
               LABELS shallowwater
               COMPILE_DEFINITIONS TYPETAG=RoughChannel LINEARSOLVER=AMGBiCGSTABIstlSolver
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_roughchannel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/roughchannel-00001.vtu
                             --zeroThreshold {"velocityY":1e-14}
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_roughchannel params.input
                             -Problem.Name roughchannel")

# convergence towards the normal flow
dumux_add_test(NAME test_shallowwater_roughchannel_convergence
               TARGET test_shallowwater_roughchannel
               LABELS shallowwater
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_shallowwater_roughchannel)

# with the second-order flux the normal flow is a discrete steady state (linear bed and free surface,
# reconstructed boundary states), so the errors are at round-off level on all grids; AMG-preconditioned
# BiCGSTAB does not converge for its Jacobian at time steps above about one second
dumux_add_test(NAME test_shallowwater_roughchannel_secondorder
               SOURCES main.cc
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=RoughChannelSecondOrder LINEARSOLVER=UMFPackIstlSolver
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_roughchannel_secondorder
               CMD_ARGS params.input -Problem.Name roughchannel_secondorder -Problem.ErrorTolerance 1e-12)

dumux_add_test(NAME test_shallowwater_roughchannel_secondorder_coarse
               TARGET test_shallowwater_roughchannel_secondorder
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_roughchannel_secondorder
               CMD_ARGS params.input -Problem.Name roughchannel_secondorder_coarse -Grid.Cells0 50
                        -Problem.ErrorTolerance 1e-12)
