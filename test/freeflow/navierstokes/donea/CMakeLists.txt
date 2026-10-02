# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)
dune_symlink_to_source_files(FILES convergencetest.py)

dumux_add_test(NAME test_ff_stokes_donea_taylorhood
               LABELS freeflow navierstokes donea
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND AND DUMUX_HAVE_GRIDFORMAT )"
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumFE
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ2FEModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MASS_MODEL=NavierStokesMassOnePFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=PQ1FEModel
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=true
               COMMAND ./convergencetest.py
               CMD_ARGS test_ff_stokes_donea_taylorhood params_taylorhood.input)

dumux_add_test(NAME test_ff_stokes_donea_momentum
               LABELS freeflow navierstokes donea
               SOURCES main_momentum.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_fcstaggered_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum -Problem.Name test_ff_stokes_donea_momentum")

dumux_add_test(NAME test_ff_stokes_donea_momentum_parallel
               LABELS freeflow navierstokes donea
               TARGET test_ff_stokes_donea_momentum
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_parallel_fcstaggered_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum -Problem.Name test_ff_stokes_donea_momentum_parallel")

dumux_add_test(NAME test_ff_stokes_donea_nocaching
              LABELS freeflow navierstokes donea
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=0
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Problem.Name test_ff_stokes_donea_nocaching -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_nocaching
              LABELS freeflow navierstokes donea
              TARGET test_ff_stokes_donea_nocaching
              CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Grid.File grids/donut.dgf
                             -Problem.Name test_ff_stokes_donea_donut_nocaching -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_twisted_nocaching
              LABELS freeflow navierstokes donea
              TARGET test_ff_stokes_donea_nocaching
              CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_twisted_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Grid.File grids/donut_twisted.dgf
                             -Problem.Name test_ff_stokes_donea_donut_twisted_nocaching -Problem.UseNeumann false"
                             --ignore "connectivity")

dumux_add_test(NAME test_ff_stokes_donea
              LABELS freeflow navierstokes donea
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=1
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                             -Problem.Name test_ff_stokes_donea -Problem.UseNeumann false")

# The staggered grid needs an overlap (Grid.Overlap in params.input), which YaspGrid provides. The residual of the
# staggered system is large in absolute terms, so Newton checks the residual reduction.
dumux_add_test(NAME test_ff_stokes_donea_stokessolver_yasp
              LABELS freeflow navierstokes donea
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=1
              COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
              COMPILE_DEFINITIONS USE_STOKES_SOLVER=1
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_stokessolver_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_stokessolver_yasp-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_stokessolver_yasp params.input
                             -Problem.Name test_ff_stokes_donea_stokessolver_yasp -Problem.UseNeumann false
                             -LinearSolver.Type gmres -LinearSolver.GMResRestart 500 -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.EnableAbsoluteResidualCriterion false -Newton.EnableResidualCriterion true
                             -Newton.ResidualReduction 1e-6 -Newton.MaxSteps 1")

dumux_add_test(NAME test_ff_stokes_donea_stokessolver_yasp_parallel
              TARGET test_ff_stokes_donea_stokessolver_yasp
              LABELS freeflow navierstokes donea
              CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_stokessolver_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_stokessolver_yasp_parallel-00001.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_stokessolver_yasp params.input
                             -Problem.Name test_ff_stokes_donea_stokessolver_yasp_parallel -Problem.UseNeumann false
                             -LinearSolver.Type gmres -LinearSolver.GMResRestart 500 -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.EnableAbsoluteResidualCriterion false -Newton.EnableResidualCriterion true
                             -Newton.ResidualReduction 1e-6 -Newton.MaxSteps 1")

dumux_add_test(NAME test_ff_navierstokes_donea_stokessolver_yasp
              TARGET test_ff_stokes_donea_stokessolver_yasp
              LABELS freeflow navierstokes donea
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_donea_stokessolver_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_donea_stokessolver_yasp-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_stokessolver_yasp params.input
                             -Problem.Name test_ff_navierstokes_donea_stokessolver_yasp -Problem.UseNeumann false
                             -Problem.EnableInertiaTerms true
                             -LinearSolver.Type gmres -LinearSolver.GMResRestart 500 -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.EnableAbsoluteResidualCriterion false -Newton.EnableResidualCriterion true
                             -Newton.ResidualReduction 1e-6")

dumux_add_test(NAME test_ff_navierstokes_donea_stokessolver_yasp_parallel
              TARGET test_ff_stokes_donea_stokessolver_yasp
              LABELS freeflow navierstokes donea
              CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_donea_stokessolver_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_navierstokes_donea_stokessolver_yasp_parallel-00001.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_stokessolver_yasp params.input
                             -Problem.Name test_ff_navierstokes_donea_stokessolver_yasp_parallel -Problem.UseNeumann false
                             -Problem.EnableInertiaTerms true
                             -LinearSolver.Type gmres -LinearSolver.GMResRestart 500 -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.EnableAbsoluteResidualCriterion false -Newton.EnableResidualCriterion true
                             -Newton.ResidualReduction 1e-6")

dumux_add_test(NAME test_ff_navierstokes_donea_stokessolver_yasp_parallel_overlap2
              TARGET test_ff_stokes_donea_stokessolver_yasp
              LABELS freeflow navierstokes donea
              CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS       --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_donea_stokessolver_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_navierstokes_donea_stokessolver_yasp_parallel_overlap2-00001.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_stokessolver_yasp params.input
                             -Problem.Name test_ff_navierstokes_donea_stokessolver_yasp_parallel_overlap2 -Problem.UseNeumann false
                             -Problem.EnableInertiaTerms true -Grid.Overlap 2
                             -LinearSolver.Type gmres -LinearSolver.GMResRestart 500 -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.EnableAbsoluteResidualCriterion false -Newton.EnableResidualCriterion true
                             -Newton.ResidualReduction 1e-6")

dumux_add_test(NAME test_ff_stokes_donea_box_stabilized
               LABELS freeflow navierstokes donea
               SOURCES main.cc
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumBox
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_box_stabilized-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_box_stabilized-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_box_stabilized params.input
                                 -Problem.Name test_ff_stokes_donea_box_stabilized -Problem.UseNeumann true
                                 -Problem.AddBoxStabilization true")

dumux_add_test(NAME test_ff_stokes_donea_donut
               LABELS freeflow navierstokes donea
               TARGET test_ff_stokes_donea
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS       --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut-00001.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                               -Grid.File grids/donut.dgf
                               -Problem.Name test_ff_stokes_donea_donut -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_twisted
               LABELS freeflow navierstokes donea
               TARGET test_ff_stokes_donea
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS       --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_twisted-00001.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                               -Grid.File grids/donut_twisted.dgf
                               -Problem.Name test_ff_stokes_donea_donut_twisted -Problem.UseNeumann false"
                               --ignore "connectivity")

# tests with diamond scheme
dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_quad
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_quad_fcdiamond_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_quad params.input
                             -Problem.Name test_ff_stokes_donea_momentum_diamond_quad -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")


dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_quad_parallel
               TARGET test_ff_stokes_donea_momentum_diamond_quad
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_diamond_quad_parallel_fcdiamond_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_quad
                             -Problem.Name test_ff_stokes_donea_momentum_diamond_quad_parallel -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_simplex_fcdiamond_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_simplex params.input
                             -Problem.Name test_ff_stokes_donea_momentum_diamond_simplex -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_diamond_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_diamond_simplex_parallel_fcdiamond_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_simplex params.input
                             -Problem.Name test_ff_stokes_donea_momentum_diamond_simplex_parallel -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

# tests with pq1 nonconforming fe
dumux_add_test(NAME test_ff_stokes_donea_momentum_fe_pq1nonconforming_cube
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ1NonconformingFEModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumFE
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS KEEP_DEFAULT_GRIDVARIABLES=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_fe_pq1nonconforming_cube-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_fe_pq1nonconforming_cube_fcdiamond_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_fe_pq1nonconforming_cube params.input
                           -Problem.Name donea_momentum_fe_pq1nonconforming_cube -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

# tests with pq1bubble
dumux_add_test(NAME test_ff_stokes_donea_momentum_pq1bubble_quad
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ1BubbleModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_quad_pq1bubble_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_quad params.input
                             -Problem.Name test_ff_stokes_donea_momentum_pq1bubble_quad -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_pq1bubble_quad_parallel
               TARGET test_ff_stokes_donea_momentum_pq1bubble_quad
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_pq1bubble_quad_parallel_pq1bubble_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_quad params.input
                             -Problem.Name test_ff_stokes_donea_momentum_pq1bubble_quad_parallel -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_pq1bubble_simplex
                SOURCES main_momentum.cc
                LABELS freeflow navierstokes donea
                CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
                COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
                COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1Bubble
                COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
                COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
                COMMAND ${DUMUX_RUNTEST}
                CMD_ARGS    --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_simplex-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_simplex_pq1bubble_1.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_simplex params.input
                            -Problem.Name test_ff_stokes_donea_momentum_pq1bubble_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_pq1bubble_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_pq1bubble_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_pq1bubble_simplex_parallel_pq1bubble_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_simplex params.input
                             -Problem.Name test_ff_stokes_donea_momentum_pq1bubble_simplex_parallel -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_box_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_box_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_box_simplex_box_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_box_simplex params.input
                           -Problem.Name test_ff_stokes_donea_momentum_box_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_box_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_box_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy --ignore "rank"
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_box_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_box_simplex_parallel_box_1.pvtu
                           --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_box_simplex params.input
                           -Problem.Name test_ff_stokes_donea_momentum_box_simplex_parallel -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1BubbleHybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex_pq1bubble_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex_parallel_pq1bubble_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex params.input
                             -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex_parallel -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_fe_pq1bubble_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ1BubbleFEModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumFE
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS KEEP_DEFAULT_GRIDVARIABLES=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_fe_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_fe_pq1bubble_simplex_pq1bubble_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_fe_pq1bubble_simplex params.input
                           -Problem.Name donea_momentum_fe_pq1bubble_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_fe_pq2_cube
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ2FEModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumFE
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS KEEP_DEFAULT_GRIDVARIABLES=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_fe_pq2_cube-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_fe_pq2_cube_pq2_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_fe_pq2_cube params.input
                           -Problem.Name donea_momentum_fe_pq2_cube -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_quadrature_pq1bubble_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1BubbleHybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS QUADRATURE_RULE=Dumux::QuadratureRules::DuneQuadrature<2>
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_quadrature_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_quadrature_pq1bubble_simplex_pq1bubble_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_quadrature_pq1bubble_simplex params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_quadrature_pq1bubble_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

# The Stokes solver tests keep the linear tolerance above the round-off floor of GMRES, where parallel
# iteration counts become erratic; the absolute Newton tolerance matches it for a single Newton step.
dumux_add_test(NAME test_ff_stokes_donea_pq1bubble_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS USE_STOKES_SOLVER=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex params.input
                           -Problem.Name test_ff_stokes_donea_pq1bubble_simplex -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -LinearSolver.Type gmres
                           -LinearSolver.GMResRestart 500
                           -LinearSolver.MaxIterations 500
                           -LinearSolver.ResidualReduction 1e-10
                           -Newton.MaxAbsoluteResidual 1e-9
                           -Newton.MaxSteps 1")

dumux_add_test(NAME test_ff_stokes_donea_pq1bubble_simplex_parallel
               TARGET test_ff_stokes_donea_pq1bubble_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_pq1bubble_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_ff_stokes_donea_pq1bubble_simplex_parallel-00001.pvtu
                             --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex params.input
                             -Problem.Name test_ff_stokes_donea_pq1bubble_simplex_parallel
                             -Problem.UseNeumann false
                             -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                             -LinearSolver.Type gmres
                             -LinearSolver.GMResRestart 500
                             -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.MaxAbsoluteResidual 1e-9
                             -Newton.MaxSteps 1")

# on four processes several of them meet at cross points; the symmetric preconditioner applies both coupling blocks
dumux_add_test(NAME test_ff_stokes_donea_pq1bubble_simplex_parallel_symmetric
               TARGET test_ff_stokes_donea_pq1bubble_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_pq1bubble_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_pq1bubble_simplex_parallel_symmetric-00001.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex params.input
                             -Problem.Name test_ff_stokes_donea_pq1bubble_simplex_parallel_symmetric
                             -Problem.UseNeumann false
                             -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                             -LinearSolver.Type gmres
                             -LinearSolver.Preconditioner.Mode Symmetric
                             -LinearSolver.GMResRestart 500
                             -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.MaxAbsoluteResidual 1e-9
                             -Newton.MaxSteps 1")

# YaspGrid has an overlap (Grid.Overlap in params.input), so the parallel run uses an overlapping decomposition
dumux_add_test(NAME test_ff_stokes_donea_hybrid_pq1bubble_yasp
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1BubbleHybrid
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS USE_STOKES_SOLVER=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_pq1bubble_yasp-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_yasp-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_yasp params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_pq1bubble_yasp
                           -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -LinearSolver.Type gmres
                           -LinearSolver.GMResRestart 500
                           -LinearSolver.MaxIterations 500
                           -LinearSolver.ResidualReduction 1e-10
                           -Newton.MaxAbsoluteResidual 1e-9
                           -Newton.MaxSteps 1")

dumux_add_test(NAME test_ff_stokes_donea_hybrid_pq1bubble_yasp_parallel
               TARGET test_ff_stokes_donea_hybrid_pq1bubble_yasp
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_pq1bubble_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_ff_stokes_donea_hybrid_pq1bubble_yasp_parallel-00001.pvtu
                             --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_yasp params.input
                             -Problem.Name test_ff_stokes_donea_hybrid_pq1bubble_yasp_parallel
                             -Problem.UseNeumann true
                             -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                             -LinearSolver.Type gmres
                             -LinearSolver.GMResRestart 500
                             -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.MaxAbsoluteResidual 1e-9
                             -Newton.MaxSteps 1")

dumux_add_test(NAME test_ff_stokes_donea_hybrid_pq1bubble_yasp_parallel_symmetric
               TARGET test_ff_stokes_donea_hybrid_pq1bubble_yasp
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_pq1bubble_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_hybrid_pq1bubble_yasp_parallel_symmetric-00001.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_yasp params.input
                             -Problem.Name test_ff_stokes_donea_hybrid_pq1bubble_yasp_parallel_symmetric
                             -Problem.UseNeumann true
                             -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                             -LinearSolver.Type gmres
                             -LinearSolver.Preconditioner.Mode Symmetric
                             -LinearSolver.GMResRestart 500
                             -LinearSolver.MaxIterations 500
                             -LinearSolver.ResidualReduction 1e-10
                             -Newton.MaxAbsoluteResidual 1e-9
                             -Newton.MaxSteps 1")


dumux_add_test(NAME test_ff_stokes_donea_hybrid_pq1bubble_cube
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=cube
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1BubbleHybrid
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_pq1bubble_cube-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_cube-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_cube params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_pq1bubble_cube -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_hybrid_pq1bubble_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1BubbleHybrid
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_simplex-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_simplex params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_pq1bubble_simplex -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_hybrid_quadrature_pq1bubble_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1BubbleHybrid
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS QUADRATURE_RULE_MOMENTUM=Dumux::QuadratureRules::DuneQuadrature<2>
               COMPILE_DEFINITIONS QUADRATURE_RULE_MASS=Dumux::QuadratureRules::MidpointQuadrature
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_quadrature_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq1bubble_simplex-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq1bubble_simplex params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_quadrature_pq1bubble_simplex -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq2_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ2Hybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq2_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_simplex_pq2_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_simplex params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq2_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq2_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_hybrid_pq2_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq2_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_hybrid_pq2_simplex_parallel_pq2_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_simplex params.input
                             -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq2_simplex_parallel -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq2_cube
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ2Hybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq2_cube-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_cube_pq2_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_cube params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq2_cube -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

# same setup as above but with the grid variables cache disabled;
# has to reproduce the cached result exactly
dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq2_cube_nocaching
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ2Hybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS ENABLEGRIDVARSCACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq2_cube-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_cube_nocaching_pq2_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_cube_nocaching params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq2_cube_nocaching -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq2_cube_parallel
               TARGET test_ff_stokes_donea_momentum_hybrid_pq2_cube
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq2_cube-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_momentum_hybrid_pq2_cube_parallel_pq2_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_cube params.input
                             -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq2_cube_parallel -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

# PQ3Hybrid momentum-only tests (cubic Lagrange)
dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq3_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND AND DUMUX_HAVE_GRIDFORMAT )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ3Hybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq3_simplex_ho-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_simplex_ho_pq3.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_simplex params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq3_simplex
                           -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -Grid.Cells \"8 8\"")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq3_cube
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND AND DUMUX_HAVE_GRIDFORMAT )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ3Hybrid
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq3_cube_ho-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_cube_ho_pq3.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_cube params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq3_cube
                           -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -Grid.Cells \"8 8\"")

# PQ3 parallel momentum tests
dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq3_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_hybrid_pq3_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND AND DUMUX_HAVE_GRIDFORMAT )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy --ignore "process rank"
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq3_simplex_ho-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_simplex_parallel_ho_pq3.pvtu
                           --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_simplex params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq3_simplex_parallel
                           -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -Grid.Cells \"8 8\"")

dumux_add_test(NAME test_ff_stokes_donea_momentum_hybrid_pq3_cube_parallel
               TARGET test_ff_stokes_donea_momentum_hybrid_pq3_cube
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND AND DUMUX_HAVE_GRIDFORMAT )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy --ignore "process rank"
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_hybrid_pq3_cube_ho-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_cube_parallel_ho_pq3.pvtu
                           --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq3_cube params.input
                           -Problem.Name test_ff_stokes_donea_momentum_hybrid_pq3_cube_parallel
                           -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -Grid.Cells \"8 8\"")

# PQ3 fe momentum tests
dumux_add_test(NAME test_ff_stokes_donea_momentum_fe_pq3_cube
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND AND DUMUX_HAVE_GRIDFORMAT )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ3FEModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumFE
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS KEEP_DEFAULT_GRIDVARIABLES=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_fe_pq3_cube_ho-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_fe_pq3_cube_ho_pq3.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_fe_pq3_cube params.input
                           -Problem.Name test_ff_stokes_donea_momentum_fe_pq3_cube
                           -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -Grid.Cells \"8 8\"")


dumux_add_test(NAME test_ff_stokes_donea_hybrid_quadrature_pq2_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ2Hybrid
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS QUADRATURE_RULE_MOMENTUM=Dumux::QuadratureRules::DuneQuadrature<4>
               COMPILE_DEFINITIONS QUADRATURE_RULE_MASS=Dumux::QuadratureRules::DuneQuadrature<4>
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMPILE_DEFINITIONS USE_STOKES_SOLVER=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_quadrature_pq2_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq2_simplex-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq2_simplex params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_quadrature_pq2_simplex -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -LinearSolver.Type gmres -LinearSolver.SymmetrizeDirichlet false -LinearSolver.DirectSolverForVelocity true")

# p2-p1 (Taylor-Hood) Stokes solver in parallel: iterative velocity AMG instead of the
# sequential-only direct velocity solver used by the test above.
dumux_add_test(NAME test_ff_stokes_donea_hybrid_quadrature_pq2_simplex_parallel
               TARGET test_ff_stokes_donea_hybrid_quadrature_pq2_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy --ignore "process rank"
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_quadrature_pq2_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_ff_stokes_donea_hybrid_quadrature_pq2_simplex_parallel-00001.pvtu
                           --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq2_simplex params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_quadrature_pq2_simplex_parallel -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -LinearSolver.Type gmres -LinearSolver.SymmetrizeDirichlet false -LinearSolver.DirectSolverForVelocity false
                           -LinearSolver.GMResRestart 500 -LinearSolver.MaxIterations 500
                           -LinearSolver.ResidualReduction 1e-10 -Newton.MaxAbsoluteResidual 1e-9")

dumux_add_test(NAME test_ff_stokes_donea_hybrid_quadrature_pq2_simplex_parallel_symmetric
               TARGET test_ff_stokes_donea_hybrid_quadrature_pq2_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy --ignore "process rank"
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_quadrature_pq2_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/s0004-test_ff_stokes_donea_hybrid_quadrature_pq2_simplex_parallel_symmetric-00001.pvtu
                           --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq2_simplex params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_quadrature_pq2_simplex_parallel_symmetric -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -LinearSolver.Type gmres -LinearSolver.Preconditioner.Mode Symmetric -LinearSolver.SymmetrizeDirichlet false -LinearSolver.DirectSolverForVelocity false
                           -LinearSolver.GMResRestart 500 -LinearSolver.MaxIterations 500
                           -LinearSolver.ResidualReduction 1e-10 -Newton.MaxAbsoluteResidual 1e-9")

dumux_add_test(NAME test_ff_stokes_donea_hybrid_quadrature_pq2_cube
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=cube
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ2Hybrid
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS QUADRATURE_RULE_MOMENTUM=Dumux::QuadratureRules::DuneQuadrature<4>
               COMPILE_DEFINITIONS QUADRATURE_RULE_MASS=Dumux::QuadratureRules::DuneQuadrature<4>
               COMPILE_DEFINITIONS ENABLEFLUXVARSCACHING=0
               COMPILE_DEFINITIONS USE_STOKES_SOLVER=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_quadrature_pq2_cube-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq2_cube-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_quadrature_pq2_cube params.input
                           -Problem.Name test_ff_stokes_donea_hybrid_quadrature_pq2_cube -Problem.UseNeumann true
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -LinearSolver.Type gmres -LinearSolver.SymmetrizeDirichlet false -LinearSolver.DirectSolverForVelocity true")

# the coupled system solved directly with MUMPS: sequentially, and in parallel on a non-overlapping
# (ALUGrid) and an overlapping (YaspGrid) decomposition, where four processes meet at cross points
dumux_add_test(NAME test_ff_stokes_donea_pq1bubble_simplex_mumps
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( MUMPS_FOUND AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1Bubble
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS USE_MUMPS_SOLVER=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex_mumps-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex_mumps params.input
                           -Problem.Name test_ff_stokes_donea_pq1bubble_simplex_mumps -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                           -Newton.MaxAbsoluteResidual 1e-9 -Newton.MaxSteps 1")

foreach(np 2 4)
    dumux_add_test(NAME test_ff_stokes_donea_pq1bubble_simplex_mumps_parallel_${np}
                   TARGET test_ff_stokes_donea_pq1bubble_simplex_mumps
                   LABELS freeflow navierstokes donea
                   CMAKE_GUARD "( MUMPS_FOUND AND dune-alugrid_FOUND )"
                   COMMAND ${DUMUX_RUNTEST}
                   CMD_ARGS      --script fuzzy --ignore "process rank"
                                 --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_pq1bubble_simplex-reference.vtu
                                         ${CMAKE_CURRENT_BINARY_DIR}/s000${np}-test_ff_stokes_donea_pq1bubble_simplex_mumps_parallel_${np}-00001.pvtu
                                 --command "${MPIEXEC} -np ${np} ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex_mumps params.input
                                 -Problem.Name test_ff_stokes_donea_pq1bubble_simplex_mumps_parallel_${np} -Problem.UseNeumann false
                                 -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                                 -Newton.MaxAbsoluteResidual 1e-9 -Newton.MaxSteps 1")
endforeach()

dumux_add_test(NAME test_ff_stokes_donea_hybrid_pq1bubble_yasp_mumps_parallel
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD MUMPS_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::YaspGrid<2>
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=DoneaTestMomentumPQ1BubbleHybrid
               COMPILE_DEFINITIONS TYPETAG_MASS=DoneaTestMassBox
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS USE_MUMPS_SOLVER=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "process rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_hybrid_pq1bubble_yasp-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_ff_stokes_donea_hybrid_pq1bubble_yasp_mumps_parallel-00001.pvtu
                             --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_hybrid_pq1bubble_yasp_mumps_parallel params.input
                             -Problem.Name test_ff_stokes_donea_hybrid_pq1bubble_yasp_mumps_parallel -Problem.UseNeumann true
                             -FreeFlow.EnableUnsymmetrizedVelocityGradient false
                             -Newton.MaxAbsoluteResidual 1e-9 -Newton.MaxSteps 1")
