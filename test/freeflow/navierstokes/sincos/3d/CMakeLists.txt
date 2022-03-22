dune_symlink_to_source_files(FILES grids)
add_input_file_links()

dumux_add_test(NAME test_ff_stokes_trigonometric_momentum_diamond_cube
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes test_ff_stokes_trigonometric_momentum_diamond_cube
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_trigonometric_momentum_diamond_cube
               CMD_ARGS params.input)

dumux_add_test(NAME test_ff_stokes_trigonometric_momentum_diamond_simplex
               SOURCES main_momentum.cc
               LABELS freeflow navierstokes test_ff_stokes_trigonometric_momentum_diamond_simplex
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_trigonometric_momentum_diamond_simplex
               CMD_ARGS params.input)

# full Stokes equations
dumux_add_test(NAME test_ff_stokes_trigonometric_diamond_tpfa_cube
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=CCTpfaModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_trigonometric_diamond_tpfa_cube
               CMD_ARGS params.input)

dumux_add_test(NAME test_ff_stokes_trigonometric_diamond_tpfa_simplex
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=CCTpfaModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_trigonometric_diamond_tpfa_simplex
               CMD_ARGS params.input)

dumux_add_test(NAME test_ff_stokes_trigonometric_diamond_box_cube
               SOURCES main.cc
               LABELS freeflow navierstokes donea
               CMAKE_GUARD dune-alugrid_FOUND
               COMPILE_DEFINITIONS GRIDTYPE=Dune::ALUGrid<3,3,Dune::cube,Dune::nonconforming>
               COMPILE_DEFINITIONS MOMENTUM_DISCRETIZATION_MODEL=FaceCenteredDiamondModel
               COMPILE_DEFINITIONS MASS_DISCRETIZATION_MODEL=BoxModel
               COMPILE_DEFINITIONS NAVIER_STOKES_MODEL=NavierStokesMomentumDiamond
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_trigonometric_diamond_box_cube
               CMD_ARGS params.input)