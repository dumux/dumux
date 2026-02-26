# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_ff_stokes_donea_momentum
               LABELS freeflow navierstokes donea
               SOURCES main_momentum.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_fcstaggered_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum -Problem.Name donea_momentum")

dumux_add_test(NAME test_ff_stokes_donea_momentum_parallel
               LABELS freeflow navierstokes donea
               TARGET test_ff_stokes_donea_momentum
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-donea_momentum_fcstaggered_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum -Problem.Name donea_momentum")

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

dumux_add_test(NAME test_ff_stokes_donea_box_stabilized
               LABELS freeflow navierstokes donea
               SOURCES main.cc
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
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
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_quad_fcdiamond_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_quad params.input
                             -Problem.Name donea_momentum_quad -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_quad_parallel
               TARGET test_ff_stokes_donea_momentum_diamond_quad
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-donea_momentum_quad_fcdiamond_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_quad
                             -Problem.Name donea_momentum_quad -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_fcdiamond_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_simplex params.input
                             -Problem.Name donea_momentum_simplex -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_diamond_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-donea_momentum_simplex_fcdiamond_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_simplex params.input
                             -Problem.Name donea_momentum_simplex -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

# tests with pq1bubble
dumux_add_test(NAME test_ff_stokes_donea_momentum_pq1bubble_quad
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ1BubbleModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_quad_pq1bubble_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_quad params.input
                             -Problem.Name donea_momentum_quad -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_pq1bubble_simplex_old_interface
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ1BubbleModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_old_interface_pq1bubble_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_simplex_old_interface params.input
                             -Problem.Name donea_momentum_simplex_old_interface -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

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
                                    ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_pq1bubble_1.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_simplex params.input
                            -Problem.Name donea_momentum_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

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
                                   ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_hybrid_pq1bubble_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq1bubble_simplex params.input
                           -Problem.Name donea_momentum_simplex_hybrid -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

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
                                   ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_hybrid_quadrature_pq1bubble_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_quadrature_pq1bubble_simplex params.input
                           -Problem.Name donea_momentum_simplex_hybrid_quadrature -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_pq1bubble_simplex_old_interface
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ1BubbleModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=0
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=0
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex_old_interface-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex_old_interface params.input
                           -Problem.Name test_ff_stokes_donea_pq1bubble_simplex_old_interface -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_pq1bubble_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ1BubbleModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS    --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_pq1bubble_simplex-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex-00001.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_pq1bubble_simplex params.input
                           -Problem.Name test_ff_stokes_donea_pq1bubble_simplex -Problem.UseNeumann false
                           -FreeFlow.EnableUnsymmetrizedVelocityGradient false")


dumux_add_test(NAME test_ff_stokes_donea_hybrid_pq1bubble_cube
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=cube
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ1BubbleHybridModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
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
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ1BubbleHybridModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
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
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ1BubbleHybridModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS PQ1BUBBLE_QUADRATURE_RULE_MOMENTUM=Dumux::QuadratureRules::DuneQuadrature<2>
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
                                   ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_hybrid_pq2_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_simplex params.input
                           -Problem.Name donea_momentum_simplex_hybrid -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

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
                                   ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_cube_hybrid_pq2_1.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_hybrid_pq2_cube params.input
                           -Problem.Name donea_momentum_cube_hybrid -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_hybrid_quadrature_pq2_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=simplex
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ2HybridModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS PQ2_QUADRATURE_RULE_MOMENTUM=Dumux::QuadratureRules::DuneQuadrature<4>
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

dumux_add_test(NAME test_ff_stokes_donea_hybrid_quadrature_pq2_cube
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS ALUGRID_CELL_TYPE=cube
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=PQ2HybridModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MOMENTUM_MODEL=NavierStokesMomentumCVFE
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NEW_PROBLEM_INTERFACE=1
               COMPILE_DEFINITIONS NEW_VARIABLES_INTERFACE=1
               COMPILE_DEFINITIONS ENABLECACHING=1
               COMPILE_DEFINITIONS PQ2_QUADRATURE_RULE_MOMENTUM=Dumux::QuadratureRules::DuneQuadrature<4>
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
