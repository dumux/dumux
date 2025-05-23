# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "convergencetest.py")
dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_ff_navierstokes_sincos
              SOURCES main.cc
              LABELS freeflow
              TIMEOUT 3600
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS LINEARSOLVER=UMFPackIstlSolver<SeqLinearSolverTraits,LinearAlgebraTraitsFromAssembler<Assembler>>
              COMMAND ./convergencetest.py
              CMD_ARGS test_ff_navierstokes_sincos params.input
              -Grid.UpperRight "6.28 6.28"
              -Grid.Cells "150 150"
              -Problem.Name test_ff_navierstokes_sincos_stationary_convergence
              -Problem.PrintConvergenceTestFile true
              -Problem.IsStationary true
              -Component.LiquidKinematicViscosity 0.001)

dumux_add_test(NAME test_ff_navierstokes_sincos_instationary
              TARGET test_ff_navierstokes_sincos
              LABELS freeflow
              TIMEOUT 3600
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS LINEARSOLVER=UMFPackIstlSolver<SeqLinearSolverTraits,LinearAlgebraTraitsFromAssembler<Assembler>>
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_instationary-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_instationary-00017.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos params.input
                             -Grid.UpperRight '1 1'
                             -Grid.Cells '50 50'
                             -Problem.Name test_ff_navierstokes_sincos_instationary
                             -Problem.IsStationary false
                             -Component.LiquidKinematicViscosity 0.1")

dumux_add_test(NAME test_ff_navierstokes_sincos_uzawapreconditioner_factory
              SOURCES main.cc
              LABELS freeflow
              TIMEOUT 5000
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS LINEARSOLVER=IstlSolverFactoryBackend<LinearSolverTraits<MassGridGeometry>,LinearAlgebraTraitsFromAssembler<Assembler>>
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_instationary-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_uzawapreconditioner-00017.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_uzawapreconditioner_factory params.input
                             -Grid.UpperRight '1 1'
                             -Grid.Cells '50 50'
                             -Problem.Name test_ff_navierstokes_sincos_uzawapreconditioner
                             -Problem.IsStationary false
                             -Component.LiquidKinematicViscosity 0.1")
