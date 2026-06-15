# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input cylinder_channel.msh cylinder_channel_quad.msh
                                   dfg_benchmark2_driver.py)

# Navier-Stokes version of the test (Re=20)
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumDiamond
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassTpfa
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_stationary_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa -Problem.EnableInertiaTerms true")

# Stokes version of the test (no inertia)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_diamond_tpfa
               TARGET test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               TIMEOUT 600
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_diamond_tpfa-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_diamond_tpfa params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_diamond_tpfa -Problem.EnableInertiaTerms false
                             -LinearSolver.UseIterativeSolver true")

# Navier-Stokes version of the test (Re=20)
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassBox
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box -Problem.EnableInertiaTerms true")

# Stokes version of the test (no inertia)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box
               TARGET test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               TIMEOUT 600
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box -Problem.EnableInertiaTerms false
                             -LinearSolver.UseIterativeSolver true")

# Navier-Stokes version of the test (Re=20)
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassDiamond
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond -Problem.EnableInertiaTerms true
                             -Grid.File cylinder_channel_quad.msh")

# Stokes version of the test (no inertia)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond
               TARGET test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_diamond params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond -Problem.EnableInertiaTerms false
                             -Grid.File cylinder_channel_quad.msh")

# Stokes (Navier-Stokes without inertia) using Trilinos Amesos2 direct solver
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond_trilinos
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassDiamond
               COMPILE_DEFINITIONS GRID_TYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               LABELS freeflow navierstokes trilinos
               CMAKE_GUARD "( Trilinos_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond_trilinos
               CMD_ARGS      params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond_trilinos
                             -Problem.EnableInertiaTerms false
                             -Grid.File cylinder_channel_quad.msh
                             -LinearSolver.UseTrilinos true)

# Parallel Stokes using Trilinos Amesos2 (2 MPI ranks, PQ1Bubble + FCDiamond)
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond_trilinos_parallel
               TARGET test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond_trilinos
               LABELS freeflow navierstokes trilinos
               CMAKE_GUARD "( Trilinos_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${MPIEXEC_EXECUTABLE}
               CMD_ARGS      ${MPIEXEC_NUMPROC_FLAG} 2
                             ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond_trilinos
                             params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_diamond_trilinos_parallel
                             -Problem.EnableInertiaTerms false
                             -Grid.File cylinder_channel_quad.msh
                             -LinearSolver.UseTrilinos true)

# Stokes (Navier-Stokes without inertia) using Trilinos Amesos2, PQ1Bubble (Mini) + Box on simplices
# Mini element on triangles is inf-sup stable; this also tests the NS nonlinear path in parallel.
dumux_add_test(NAME test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassBox
               COMPILE_DEFINITIONS GRID_TYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               LABELS freeflow navierstokes trilinos
               CMAKE_GUARD "( Trilinos_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos
               CMD_ARGS      params.input
                             -Problem.Name test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos
                             -Problem.EnableInertiaTerms false
                             -LinearSolver.UseTrilinos true)

# Parallel NS using Trilinos Amesos2 (2 MPI ranks, Mini element + Box on simplices)
# Mini element is inf-sup stable on triangles; tests nonlinear NS parallel direct solve.
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box_trilinos_parallel
               TARGET test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos
               LABELS freeflow navierstokes trilinos
               CMAKE_GUARD "( Trilinos_FOUND AND dune-alugrid_FOUND )"
               COMMAND ${MPIEXEC_EXECUTABLE}
               CMD_ARGS      ${MPIEXEC_NUMPROC_FLAG} 2
                             ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos
                             params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_stationary_pq1bubble_box_trilinos_parallel
                             -Problem.EnableInertiaTerms true
                             -LinearSolver.UseTrilinos true)

# Instationary (time-dependent) Navier-Stokes (DFG cylinder, Re~100) solved with the
# matrix-free SIMPLE block preconditioner (PQ1Bubble/Mini + Box on simplices, ALUGrid).
# The time-derivative term makes the velocity operator diagonally dominant, so diag(A)^-1 is
# an accurate Schur approximation and SIMPLE is effective. This is the sequential reference run.
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_instationary_simple
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DFGChannelTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DFGChannelTestMassBox
               COMPILE_DEFINITIONS GRID_TYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               TIMEOUT 1200
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_instationary_simple-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_instationary_simple-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_instationary_simple params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_instationary_simple
                             -Problem.EnableInertiaTerms true -Problem.Instationary true
                             -Component.LiquidDynamicViscosity 2e-4
                             -TimeLoop.DtInitial 0.002 -TimeLoop.TEnd 0.004 -TimeLoop.MaxTimeStepSize 0.002
                             -LinearSolver.UseSimpleSolver true -LinearSolver.Type fgmres
                             -LinearSolver.MaxIterations 500 -LinearSolver.ResidualReduction 1e-7
                             -LinearSolver.Preconditioner.SchurReduction 1e-2 -LinearSolver.Preconditioner.SchurIterations 30
                             -LinearSolver.Preconditioner.PressureRelaxation 1.0
                             -Newton.MaxSteps 15")

# Parallel (2 MPI ranks) instationary NS with the matrix-free SIMPLE preconditioner.
# Compares against the same sequential reference (must reproduce the direct/sequential solve).
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_instationary_simple_parallel
               TARGET test_ff_navierstokes_dfg_benchmark_instationary_simple
               LABELS freeflow navierstokes parallel
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               TIMEOUT 1200
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_instationary_simple-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_ff_navierstokes_dfg_benchmark_instationary_simple_parallel-00002.pvtu
                             --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_dfg_benchmark_instationary_simple params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_instationary_simple_parallel
                             -Problem.EnableInertiaTerms true -Problem.Instationary true
                             -Component.LiquidDynamicViscosity 2e-4
                             -TimeLoop.DtInitial 0.002 -TimeLoop.TEnd 0.004 -TimeLoop.MaxTimeStepSize 0.002
                             -LinearSolver.UseSimpleSolver true -LinearSolver.Type fgmres
                             -LinearSolver.MaxIterations 500 -LinearSolver.ResidualReduction 1e-7
                             -LinearSolver.Preconditioner.SchurReduction 1e-2 -LinearSolver.Preconditioner.SchurIterations 30
                             -LinearSolver.Preconditioner.PressureRelaxation 1.0
                             -Newton.MaxSteps 15")

# Cross-check: the SAME instationary NS setup (same discretization, Re~100, same dt) solved in
# parallel with the Trilinos Amesos2 DIRECT solver. Compared against the same reference, so the
# matrix-free SIMPLE solve and the parallel direct solve are verified to agree.
dumux_add_test(NAME test_ff_navierstokes_dfg_benchmark_instationary_trilinos_parallel
               TARGET test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos
               LABELS freeflow navierstokes trilinos parallel
               CMAKE_GUARD "( Trilinos_FOUND AND dune-alugrid_FOUND )"
               TIMEOUT 1200
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_dfg_benchmark_instationary_simple-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_ff_navierstokes_dfg_benchmark_instationary_trilinos_parallel-00002.pvtu
                             --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_dfg_benchmark_stationary_pq1bubble_box_trilinos params.input
                             -Problem.Name test_ff_navierstokes_dfg_benchmark_instationary_trilinos_parallel
                             -Problem.EnableInertiaTerms true -Problem.Instationary true
                             -Component.LiquidDynamicViscosity 2e-4
                             -TimeLoop.DtInitial 0.002 -TimeLoop.TEnd 0.004 -TimeLoop.MaxTimeStepSize 0.002
                             -LinearSolver.UseTrilinos true
                             -Newton.MaxSteps 15")
