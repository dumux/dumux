# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_ff_stokes_donea_momentum
               LABELS freeflow navierstokes donea
               SOURCES main_momentum.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_fcstaggered_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum -Problem.Name donea_momentum")

dumux_add_test(NAME test_ff_stokes_donea_momentum_parallel
               LABELS freeflow navierstokes donea
               TARGET test_ff_stokes_donea_momentum
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy --ignore "rank"
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-donea_momentum_fcstaggered_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum -Problem.Name donea_momentum")

dumux_add_test(NAME test_ff_stokes_donea_nocaching
              LABELS freeflow navierstokes donea
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Problem.Name test_ff_stokes_donea_nocaching -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_nocaching
              LABELS freeflow navierstokes donea
              TARGET test_ff_stokes_donea_nocaching
              CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_quad_fcdiamond_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_quad params.input
                             -Problem.Name donea_momentum_quad -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_quad_parallel
               TARGET test_ff_stokes_donea_momentum_diamond_quad
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_diamond_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_fcdiamond_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_simplex params.input
                             -Problem.Name donea_momentum_simplex -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient true")

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_simplex_parallel
               TARGET test_ff_stokes_donea_momentum_diamond_simplex
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND dune-alugrid_FOUND )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_quad-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_quad_pq1bubble_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_quad params.input
                             -Problem.Name donea_momentum_quad -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient false")

dumux_add_test(NAME test_ff_stokes_donea_momentum_pq1bubble_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-alugrid_FOUND )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=PQ1BubbleModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumCVFE
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum_pq1bubble_simplex-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_simplex_pq1bubble_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_pq1bubble_simplex params.input
                             -Problem.Name donea_momentum_simplex -Problem.UseNeumann false -FreeFlow.EnableUnsymmetrizedVelocityGradient false")
