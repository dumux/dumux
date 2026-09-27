# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_input_file_links()
dumux_add_test(NAME test_mortar_darcy_1p
              LABELS multidomain multidomain_mortar
              SOURCES main.cc
              CMAKE_GUARD dune-foamgrid_FOUND)

# the mortar is imposed per trace vertex and its flux read back per trace cell, which is not
# the adjoint of the imposition, so the interface operator needs a nonsymmetric Krylov method
dumux_add_test(NAME test_mortar_darcy_1p_box
              LABELS multidomain multidomain_mortar
              SOURCES main.cc
              COMPILE_DEFINITIONS MORTAR_BOX=1
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p_box>
              CMD_ARGS -Mortar.Solver gmres)

dumux_add_test(NAME test_mortar_darcy_1p_jump
              LABELS multidomain multidomain_mortar
              SOURCES main.cc
              COMPILE_DEFINITIONS PERMEABILITY_JUMP=1
              CMAKE_GUARD dune-foamgrid_FOUND)

dumux_add_test(NAME test_mortar_darcy_1p_2x1
              LABELS multidomain multidomain_mortar
              SOURCES main.cc
              COMPILE_DEFINITIONS SUBDOMAINS_2X1=1
              CMAKE_GUARD dune-foamgrid_FOUND)

dumux_add_test(NAME test_mortar_darcy_1p_2x1_prec
              LABELS multidomain multidomain_mortar
              TARGET test_mortar_darcy_1p_2x1
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p_2x1>
              CMD_ARGS -Mortar.Preconditioner interface -Mortar.CheckPreconditioner true)

dumux_add_test(NAME test_mortar_darcy_1p_jump_prec
              LABELS multidomain multidomain_mortar
              TARGET test_mortar_darcy_1p_jump
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p_jump>
              CMD_ARGS -Mortar.Preconditioner interface -Mortar.CheckPreconditioner true)

dumux_add_test(NAME test_mortar_darcy_1p_prec
              LABELS multidomain multidomain_mortar
              TARGET test_mortar_darcy_1p
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p>
              CMD_ARGS -Mortar.Preconditioner interface -Mortar.CheckPreconditioner true)

dumux_add_test(NAME test_mortar_darcy_1p_2x1_flux
              LABELS multidomain multidomain_mortar
              TARGET test_mortar_darcy_1p_2x1
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p_2x1>
              CMD_ARGS -Mortar.Variant flux)

dumux_add_test(NAME test_mortar_darcy_1p_2x1_flux_prec
              LABELS multidomain multidomain_mortar
              TARGET test_mortar_darcy_1p_2x1
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p_2x1>
              CMD_ARGS -Mortar.Variant flux -Mortar.Preconditioner interface -Mortar.CheckPreconditioner true)

dumux_add_test(NAME test_mortar_darcy_1p_jump_flux_prec
              LABELS multidomain multidomain_mortar
              TARGET test_mortar_darcy_1p_jump
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p_jump>
              CMD_ARGS -Mortar.Variant flux -Mortar.Preconditioner interface -Mortar.CheckPreconditioner true)

dumux_add_test(NAME test_mortar_darcy_1p_flux_prec
              LABELS multidomain multidomain_mortar
              TARGET test_mortar_darcy_1p
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND $<TARGET_FILE:test_mortar_darcy_1p>
              CMD_ARGS -Mortar.Variant flux -Mortar.Preconditioner interface -Mortar.CheckPreconditioner true
                       -Mortar.ResidualReduction 1e-9 -LinearSolver.ResidualReduction 1e-13)
