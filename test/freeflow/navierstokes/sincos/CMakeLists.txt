dune_symlink_to_source_files(FILES "convergencetest.py")
dune_symlink_to_source_files(FILES "params.input" "params_momentum.input")
dune_symlink_to_source_files(FILES "run_simulations.py")

add_subdirectory(3d)

dumux_add_test(NAME test_ff_navierstokes_sincos
              SOURCES main.cc
              LABELS freeflow
              TIMEOUT 3600
              CMAKE_GUARD HAVE_UMFPACK
              COMPILE_DEFINITIONS LINEARSOLVER=UMFPackBackend
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
              COMPILE_DEFINITIONS LINEARSOLVER=UMFPackBackend
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
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
              CMAKE_GUARD "( HAVE_UMFPACK AND ( ( DUNE_ISTL_VERSION VERSION_GREATER 2.7 ) OR ( DUNE_ISTL_VERSION VERSION_EQUAL 2.7 ) ) )"
              COMPILE_DEFINITIONS LINEARSOLVER=IstlSolverFactoryBackend<LinearSolverTraits<MassGridGeometry>>
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_sincos_instationary-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_uzawapreconditioner-00017.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_uzawapreconditioner_factory params.input
                             -Grid.UpperRight '1 1'
                             -Grid.Cells '50 50'
                             -Problem.Name test_ff_navierstokes_sincos_uzawapreconditioner
                             -Problem.IsStationary false
                             -Component.LiquidKinematicViscosity 0.1")

# tests with diamond scheme
# navier-stokes momentum only
dumux_add_test(NAME test_ff_navierstokes_sincos_momentum_diamond_quad
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes sincos
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_momentum_diamond_quad
               CMD_ARGS params_momentum.input)

dumux_add_test(NAME test_ff_navierstokes_sincos_momentum_diamond_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes sincos
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_momentum_diamond_simplex
               CMD_ARGS params_momentum.input)

dumux_add_test(NAME test_ff_navierstokes_sincos_momentum_staggered_quad
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes sincos
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredStaggeredModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentum
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_momentum_staggered_quad
               CMD_ARGS params_momentum.input)

# full Navier-Stokes equations
dumux_add_test(NAME test_ff_navierstokes_sincos_diamond_tpfa_quad
               SOURCES main.cc
               LABELS freeflow navierstokes sincos
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS LINEARSOLVER=UMFPackBackend
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=CCTpfaModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_diamond_tpfa_quad
               CMD_ARGS params.input -Problem.UseNeumann true)

dumux_add_test(NAME test_ff_navierstokes_sincos_diamond_box_quad
               SOURCES main.cc
               LABELS freeflow navierstokes sincos
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS LINEARSOLVER=UMFPackBackend
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_diamond_box_quad
               CMD_ARGS params.input)

dumux_add_test(NAME test_ff_navierstokes_sincos_diamond_tpfa_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes sincos
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS LINEARSOLVER=UMFPackBackend
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=CCTpfaModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_diamond_tpfa_simplex
               CMD_ARGS params.input)

dumux_add_test(NAME test_ff_navierstokes_sincos_diamond_box_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes sincos
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS LINEARSOLVER=UMFPackBackend
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_sincos_diamond_box_simplex
               CMD_ARGS params.input)
