add_input_file_links()
dune_symlink_to_source_files(FILES grids)

# tests with staggered
dumux_add_test(NAME test_ff_stokes_donea_momentum
               LABELS freeflow navierstokes donea
               SOURCES main_momentum.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/donea_momentum_fcstaggered_1.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum")

dumux_add_test(NAME test_ff_stokes_donea_momentum_parallel
               LABELS freeflow navierstokes donea
               TARGET test_ff_stokes_donea_momentum
               CMAKE_GUARD "( HAVE_UMFPACK AND HAVE_MPI AND ( DUNE_ISTL_VERSION VERSION_GREATER_EQUAL 2.8 ) )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy --zeroThreshold {"rank":100}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_momentum-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/s0004-donea_momentum_fcstaggered_1.pvtu
                             --command "${MPIEXEC} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum")

dumux_add_test(NAME test_ff_stokes_donea_nocaching
              LABELS freeflow navierstokes donea
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=false
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching_fcstaggered-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Problem.Name test_ff_stokes_donea_nocaching -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut_nocaching
              LABELS freeflow navierstokes donea
              TARGET test_ff_stokes_donea_nocaching
              CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_nocaching_fcstaggered-00001.vtu
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
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_twisted_nocaching_fcstaggered-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_nocaching params.input
                             -Grid.File grids/donut_twisted.dgf
                             -Problem.Name test_ff_stokes_donea_donut_twisted_nocaching -Problem.UseNeumann false"
                             --zeroThreshold {"connectivity":1e6})

dumux_add_test(NAME test_ff_stokes_donea
              LABELS freeflow navierstokes donea
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS ENABLECACHING=true
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_fcstaggered-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                             -Problem.Name test_ff_stokes_donea -Problem.UseNeumann false")

dumux_add_test(NAME test_ff_stokes_donea_donut
               LABELS freeflow navierstokes donea
               TARGET test_ff_stokes_donea
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_donea_donut-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_fcstaggered-00001.vtu
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
                                       ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_donut_twisted_fcstaggered-00001.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea params.input
                               -Grid.File grids/donut_twisted.dgf
                               -Problem.Name test_ff_stokes_donea_donut_twisted -Problem.UseNeumann false"
                               --zeroThreshold {"connectivity":1e6})

# tests with diamond scheme
# Stokes momentum only
dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_quad
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_quad
               CMD_ARGS params.input -Problem.UseNeumann true)

dumux_add_test(NAME test_ff_stokes_donea_momentum_diamond_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_momentum_diamond_simplex
               CMD_ARGS params.input -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient false)

# full Stokes equations
dumux_add_test(NAME test_ff_stokes_donea_diamond_tpfa_quad
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=CCTpfaModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_diamond_tpfa_quad
               CMD_ARGS params.input -Problem.UseNeumann true)

dumux_add_test(NAME test_ff_stokes_donea_diamond_tpfa_quad_unstructured
               TARGET test_ff_stokes_donea_diamond_tpfa_quad
               LABELS freeflow navierstokes donea
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_diamond_tpfa_quad
               CMD_ARGS params.input -Problem.UseNeumann false -Grid.File ./grids/unstructured_quad.msh -Grid.Refinement 2)

dumux_add_test(NAME test_ff_stokes_donea_diamond_box_quad
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_diamond_box_quad
               CMD_ARGS params.input -Problem.UseNeumann true)

dumux_add_test(NAME test_ff_stokes_donea_diamond_tpfa_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=CCTpfaModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_diamond_tpfa_simplex
               CMD_ARGS params.input -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient false)

dumux_add_test(NAME test_ff_stokes_donea_diamond_box_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_donea_diamond_box_simplex
               CMD_ARGS params.input -Problem.UseNeumann true -FreeFlow.EnableUnsymmetrizedVelocityGradient false)
