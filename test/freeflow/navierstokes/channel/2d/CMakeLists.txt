# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

add_executable(test_ff_channel EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_ff_stokes_channel_outflow
              TARGET test_ff_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_outflow-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input -Problem.OutletCondition Outflow
                             -Problem.Name test_ff_stokes_channel_outflow")

dumux_add_test(NAME test_ff_stokes_channel_outflow_gravity
              TARGET test_ff_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --relative 1e-5
                             --ignore p
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_outflow_gravity-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input -Problem.OutletCondition Outflow
                             -Problem.EnableGravity true -Problem.Name test_ff_stokes_channel_outflow_gravity")

dumux_add_test(NAME test_ff_stokes_channel_neumann_x_dirichlet_y
              TARGET test_ff_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_nxdy-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input -Problem.OutletCondition NeumannX_DirichletY
                             -Problem.Name test_ff_stokes_channel_nxdy")

dumux_add_test(NAME test_ff_stokes_channel_neumann_x_neumann_y
              TARGET test_ff_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_neumann-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_nxny-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input -Grid.UpperRight \"2 1\" -Grid.Cells \"50 25\"
                         -Problem.Name test_ff_stokes_channel_nxny -Problem.OutletCondition NeumannX_NeumannY
                         -Problem.UseVelocityProfile true -Problem.OutletPressure 0 -Problem.IsStationary true")

dumux_add_test(NAME test_ff_stokes_channel_unconstrained_outflow
              TARGET test_ff_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_unconstrained_outflow-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_unconstrained_outflow-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input -Problem.OutletCondition UnconstrainedOutflow -Grid.UpperRight \"2 1\" -Grid.Cells \"50 25\"
                         -Problem.Name test_ff_stokes_channel_unconstrained_outflow -Problem.OutletCondition UnconstrainedOutflow -Problem.UseVelocityProfile true -Problem.OutletPressure 0")

dumux_add_test(NAME test_ff_navierstokes_channel
              TARGET test_ff_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_channel-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params_navierstokes.input
                             -Vtk.WriteFaceData 1
                             -Problem.Name test_ff_navierstokes_channel")

dumux_add_test(NAME test_ff_stokes_channel_iterative
              TARGET test_ff_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 600
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_iterative-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input -Problem.OutletCondition Outflow
                             -LinearSolver.UseIterativeSolver true -Problem.Name test_ff_stokes_channel_iterative")

# A small time step, where the storage term dominates the velocity block: with the transient pressure operator
# the Stokes solver needs fewer than 100 iterations per solve, with the default pressure mass matrix about 400,
# so the iteration limit fails if the pressure operator is not used.
dumux_add_test(NAME test_ff_stokes_channel_transient_stokessolver
              SOURCES main_stokessolver.cc
              LABELS freeflow navierstokes
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_transient_stokessolver-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_transient_stokessolver-00003.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_transient_stokessolver params.input
                             -Problem.Name test_ff_stokes_channel_transient_stokessolver
                             -TimeLoop.DtInitial 1e-3 -TimeLoop.TEnd 3e-3
                             -LinearSolver.PressureOperator Transient -LinearSolver.Type fgmres
                             -LinearSolver.Preconditioner.Mode Triangular -LinearSolver.Preconditioner.AMGForPressure true
                             -LinearSolver.Preconditioner.SSORForVelocity true -LinearSolver.MaxIterations 150
                             -LinearSolver.CheckKeptMatrix true")

dumux_add_test(NAME test_ff_stokes_channel_transient_stokessolver_parallel
              TARGET test_ff_stokes_channel_transient_stokessolver
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_MPI
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel_transient_stokessolver-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_channel_transient_stokessolver_parallel-00003.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel_transient_stokessolver params.input
                             -Problem.Name test_ff_stokes_channel_transient_stokessolver_parallel
                             -TimeLoop.DtInitial 1e-3 -TimeLoop.TEnd 3e-3
                             -LinearSolver.PressureOperator Transient -LinearSolver.Type fgmres
                             -LinearSolver.Preconditioner.Mode Triangular -LinearSolver.Preconditioner.AMGForPressure true
                             -LinearSolver.Preconditioner.SSORForVelocity true -LinearSolver.MaxIterations 150
                             -LinearSolver.CheckKeptMatrix true")

add_executable(test_ff_stokesni_channel EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_stokesni_channel PUBLIC "NONISOTHERMAL=1")

dumux_add_test(NAME test_ff_stokesni_channel_convection
              TARGET test_ff_stokesni_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokesni_channel_convection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel_convection-00005.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel params_convection.input
                             -Problem.Name test_ff_stokesni_channel_convection")

dumux_add_test(NAME test_ff_stokesni_channel_conduction
              TARGET test_ff_stokesni_channel
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokesni_channel_conduction-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel_conduction-00004.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel params_conduction.input
                             -Problem.Name test_ff_stokesni_channel_conduction"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-20})
