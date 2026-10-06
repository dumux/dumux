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

dumux_add_test(NAME test_ff_navierstokes_sincos_hybrid_pq1bubble
              SOURCES main.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS LINEARSOLVER=UMFPackIstlSolver<SeqLinearSolverTraits,LinearAlgebraTraitsFromAssembler<Assembler>>
              COMPILE_DEFINITIONS TYPETAG_MOMENTUM=SincosTestMomentumPQ1BubbleHybrid
              COMPILE_DEFINITIONS TYPETAG_MASS=SincosTestMassBox
              COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_hybrid_pq1bubble-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq1bubble-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq1bubble params.input
                             -Grid.Cells '32 32'
                             -Problem.Name test_ff_navierstokes_sincos_hybrid_pq1bubble
                             -Problem.IsStationary true
                             -Problem.PrintErrors true")

dumux_add_test(NAME test_ff_navierstokes_sincos_hybrid_pq1bubble_instationary
              TARGET test_ff_navierstokes_sincos_hybrid_pq1bubble
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_hybrid_pq1bubble_instationary-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq1bubble_instationary-00017.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq1bubble params.input
                             -Grid.UpperRight '1 1'
                             -Grid.Cells '50 50'
                             -Problem.Name test_ff_navierstokes_sincos_hybrid_pq1bubble_instationary
                             -Problem.IsStationary false
                             -Component.LiquidKinematicViscosity 0.1
                             -Problem.PrintErrors true")

dumux_add_test(NAME test_ff_navierstokes_sincos_hybrid_pq2
              SOURCES main.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS LINEARSOLVER=UMFPackIstlSolver<SeqLinearSolverTraits,LinearAlgebraTraitsFromAssembler<Assembler>>
              COMPILE_DEFINITIONS TYPETAG_MOMENTUM=SincosTestMomentumPQ2Hybrid
              COMPILE_DEFINITIONS TYPETAG_MASS=SincosTestMassBox
              COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_hybrid_pq2-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq2-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq2 params.input
                             -Grid.Cells '32 32'
                             -Problem.Name test_ff_navierstokes_sincos_hybrid_pq2
                             -Problem.IsStationary true
                             -Problem.PrintErrors true")

dumux_add_test(NAME test_ff_navierstokes_sincos_hybrid_pq2_instationary
              TARGET test_ff_navierstokes_sincos_hybrid_pq2
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_hybrid_pq2_instationary-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq2_instationary-00017.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_hybrid_pq2 params.input
                             -Grid.UpperRight '1 1'
                             -Grid.Cells '50 50'
                             -Problem.Name test_ff_navierstokes_sincos_hybrid_pq2_instationary
                             -Problem.IsStationary false
                             -Component.LiquidKinematicViscosity 0.1
                             -Problem.PrintErrors true")

dumux_add_test(NAME test_ff_navierstokes_sincos_momentum_multistage_pq1bubble
              SOURCES main_momentum_multistage.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS TYPETAG_MOMENTUM_ONLY=SincosTestMomentumOnlyPQ1BubbleHybrid
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_momentum_multistage_pq1bubble
              CMD_ARGS params.input
                       -Grid.UpperRight "1 1"
                       -Grid.Cells "8 8"
                       -Problem.Name test_ff_navierstokes_sincos_momentum_multistage_pq1bubble
                       -Problem.IsStationary false
                       -Problem.EnableInertiaTerms false
                       -Component.LiquidKinematicViscosity 0.001
                       -TimeLoop.TEnd 1.0
                       -Newton.MaxRelativeShift 1e-12
                       -Newton.EnableAbsoluteResidualCriterion false
                       -Newton.SatisfyResidualAndShiftCriterion false
                       -Newton.Verbosity 0
                       -MultiStageTest.NumStepsCoarse 4
                       -MultiStageTest.NumTimeStepRefinements 4
                       -MultiStageTest.EquivalenceTolerance 1e-8
                       -MultiStageTest.SingleSolveTolerance 1e-5
                       -LinearPDESolver.Verbosity 0
                       -MultiStageTest.RateTolerance 0.2)

dumux_add_test(NAME test_ff_navierstokes_sincos_momentum_multistage_pq2
              SOURCES main_momentum_multistage.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS TYPETAG_MOMENTUM_ONLY=SincosTestMomentumOnlyPQ2Hybrid
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_momentum_multistage_pq2
              CMD_ARGS params.input
                       -Grid.UpperRight "1 1"
                       -Grid.Cells "8 8"
                       -Problem.Name test_ff_navierstokes_sincos_momentum_multistage_pq2
                       -Problem.IsStationary false
                       -Problem.EnableInertiaTerms false
                       -Component.LiquidKinematicViscosity 0.001
                       -TimeLoop.TEnd 1.0
                       -Newton.MaxRelativeShift 1e-12
                       -Newton.EnableAbsoluteResidualCriterion false
                       -Newton.SatisfyResidualAndShiftCriterion false
                       -Newton.Verbosity 0
                       -MultiStageTest.NumStepsCoarse 4
                       -MultiStageTest.NumTimeStepRefinements 4
                       -MultiStageTest.EquivalenceTolerance 1e-8
                       -MultiStageTest.SingleSolveTolerance 1e-5
                       -LinearPDESolver.Verbosity 0
                       -MultiStageTest.RateTolerance 0.2)

dumux_add_test(NAME test_ff_navierstokes_sincos_multistage_pq1bubble
              SOURCES main_multistage.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS TYPETAG_MOMENTUM=SincosTestMomentumPQ1BubbleHybrid
              COMPILE_DEFINITIONS TYPETAG_MASS=SincosTestMassBox
              COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_multistage_pq1bubble
              CMD_ARGS params.input
                       -Grid.UpperRight "1 1"
                       -Grid.Cells "8 8"
                       -Problem.Name test_ff_navierstokes_sincos_multistage_pq1bubble
                       -Problem.IsStationary false
                       -Problem.EnableInertiaTerms false
                       -Component.LiquidKinematicViscosity 0.001
                       -TimeLoop.TEnd 1.0
                       -Newton.MaxRelativeShift 1e-12
                       -Newton.EnableAbsoluteResidualCriterion false
                       -Newton.SatisfyResidualAndShiftCriterion false
                       -Newton.Verbosity 0
                       -MultiStageTest.NumStepsCoarse 4
                       -MultiStageTest.NumTimeStepRefinements 4
                       -MultiStageTest.EquivalenceTolerance 1e-8
                       -MultiStageTest.SingleSolveTolerance 1e-5
                       -LinearPDESolver.Verbosity 0
                       -MultiStageTest.RateTolerance 0.2)

dumux_add_test(NAME test_ff_navierstokes_sincos_multistage_pq2
              SOURCES main_multistage.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS TYPETAG_MOMENTUM=SincosTestMomentumPQ2Hybrid
              COMPILE_DEFINITIONS TYPETAG_MASS=SincosTestMassBox
              COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_multistage_pq2
              CMD_ARGS params.input
                       -Grid.UpperRight "1 1"
                       -Grid.Cells "8 8"
                       -Problem.Name test_ff_navierstokes_sincos_multistage_pq2
                       -Problem.IsStationary false
                       -Problem.EnableInertiaTerms false
                       -Component.LiquidKinematicViscosity 0.001
                       -TimeLoop.TEnd 1.0
                       -Newton.MaxRelativeShift 1e-12
                       -Newton.EnableAbsoluteResidualCriterion false
                       -Newton.SatisfyResidualAndShiftCriterion false
                       -Newton.Verbosity 0
                       -MultiStageTest.NumStepsCoarse 4
                       -MultiStageTest.NumTimeStepRefinements 4
                       -MultiStageTest.EquivalenceTolerance 1e-8
                       -MultiStageTest.SingleSolveTolerance 1e-5
                       -LinearPDESolver.Verbosity 0
                       -MultiStageTest.RateTolerance 0.2)
