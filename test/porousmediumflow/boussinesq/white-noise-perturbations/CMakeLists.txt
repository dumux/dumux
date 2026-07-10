# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

# --- vorticity/streamfunction formulation ------------------------------------

# 2D: streamfunction-vorticity (1 potential eq + 1 transport)
dumux_add_test(NAME test_boussinesq_white_noise_vorticity_box
              LABELS porousmediumflow boussinesq
              SOURCES main_vorticity.cc
              COMPILE_DEFINITIONS TYPETAG=BoussinesqOneSidedRB
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_vorticity_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_white_noise_vorticity_box)

dumux_add_test(NAME test_boussinesq_white_noise_vorticity_box_convergence
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_vorticity_box
              TIMEOUT 100000
              COMMAND ${Python3_EXECUTABLE}
              CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence.py
                       ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_vorticity_box
                       5)

# 3D: vector-potential (3 potential eqs + 1 transport)
dumux_add_test(NAME test_boussinesq_white_noise_vorticity_box_3d
              LABELS porousmediumflow boussinesq
              SOURCES main_vorticity.cc
              COMPILE_DEFINITIONS TYPETAG=BoussinesqOneSidedRB3D
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_vorticity_box_3d
              CMD_ARGS params3d.input -Problem.Name test_boussinesq_white_noise_vorticity_box_3d)

# --- pressure-based formulation (2D only) ------------------------------------

dumux_add_test(NAME test_boussinesq_white_noise_pressure_box
              LABELS porousmediumflow boussinesq
              SOURCES main_pressure.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_white_noise_pressure_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_white_noise_pressure_box)
