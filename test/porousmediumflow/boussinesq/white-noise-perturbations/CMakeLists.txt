# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

# --- vorticity/streamfunction formulation ------------------------------------

# 2D: streamfunction-vorticity (1 potential eq + 1 transport).
# ALUGrid simplex, parallel + h-adaptive (grad(C) indicator replaces the former static
# y-grading, see [Adaptive] in params.input) -- the machinery from
# ../common/adaptive/, ILU+GMRES parallel solver (AMG rejected, ../NOTES.md).
dumux_add_test(NAME test_boussinesq_white_noise_vorticity_box
              LABELS porousmediumflow boussinesq
              SOURCES main_vorticity.cc
              COMPILE_DEFINITIONS TYPETAG=BoussinesqOneSidedRB
              CMAKE_GUARD dune-alugrid_FOUND
              MPI_RANKS 1 2 4
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_vorticity_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_white_noise_vorticity_box)

# Same run with dynamic load balancing switched on (>1 rank only, see ../NOTES.md):
# rebalancing only moves data between ranks and must not change the physics.
dumux_add_test(NAME test_boussinesq_white_noise_vorticity_box_rebalance
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_vorticity_box
              CMAKE_GUARD dune-alugrid_FOUND
              MPI_RANKS 2 4
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_vorticity_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_white_noise_vorticity_box_rebalance
                       -Adaptive.EnableDynamicLoadBalancing true
                       -Adaptive.RebalanceEvery 3)

# 3D: vector-potential (3 potential eqs + 1 transport).
# SERIAL_YASP_3D keeps the unchanged serial YaspGrid + UMFPack path in main_vorticity.cc
# (the 3D conforming-bisection ALUGrid variant is untested for the rebalance machinery
# and deliberately out of scope, see the file header there).
dumux_add_test(NAME test_boussinesq_white_noise_vorticity_box_3d
              LABELS porousmediumflow boussinesq
              SOURCES main_vorticity.cc
              COMPILE_DEFINITIONS TYPETAG=BoussinesqOneSidedRB3D SERIAL_YASP_3D=1
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_vorticity_box_3d
              CMD_ARGS params3d.input -Problem.Name test_boussinesq_white_noise_vorticity_box_3d)

# --- pressure-based formulation (2D only) ------------------------------------

# ALUGrid simplex, parallel + h-adaptive + optional dynamic load balancing -- same
# machinery/rationale as the vorticity target above: box/CVFE discretization is
# grid-agnostic, so common/adaptive/boxrebalance.hh and griddatatransfer.hh attach
# here unchanged.
dumux_add_test(NAME test_boussinesq_white_noise_pressure_box
              LABELS porousmediumflow boussinesq
              SOURCES main_pressure.cc
              CMAKE_GUARD dune-alugrid_FOUND
              MPI_RANKS 1 2 4
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_pressure_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_white_noise_pressure_box)

# Same run with dynamic load balancing switched on (>1 rank only, see ../NOTES.md):
# rebalancing only moves data between ranks and must not change the physics.
dumux_add_test(NAME test_boussinesq_white_noise_pressure_box_rebalance
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_pressure_box
              CMAKE_GUARD dune-alugrid_FOUND
              MPI_RANKS 2 4
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_pressure_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_white_noise_pressure_box_rebalance
                       -Adaptive.EnableDynamicLoadBalancing true
                       -Adaptive.RebalanceEvery 3)

# AMG+GMRES variant (USE_AMG_SOLVER switches the LinearSolver alias in main_pressure.cc,
# everything else identical) -- ILU+GMRES(30) hits a hard convergence wall on this operator
# at Ra=8000 (confirmed: Newton solver never recovers, halves dt down to the noise floor and
# aborts). Same operator/fix as test_boussinesq_cosine_pressure_box_amg, just a harder Ra;
# see the LinearSolver comment in main_pressure.cc.
dumux_add_test(NAME test_boussinesq_white_noise_pressure_box_amg
              LABELS porousmediumflow boussinesq
              SOURCES main_pressure.cc
              COMPILE_DEFINITIONS USE_AMG_SOLVER=1
              CMAKE_GUARD dune-alugrid_FOUND
              MPI_RANKS 1
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_pressure_box_amg
              CMD_ARGS params.input -Problem.Name test_boussinesq_white_noise_pressure_box_amg)
