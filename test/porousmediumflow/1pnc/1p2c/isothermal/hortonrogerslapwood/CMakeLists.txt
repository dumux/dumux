# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input run_benchmark.py)

# growing convection cell at Ra = 8 pi^2
dumux_add_test(NAME test_1p2c_hortonrogerslapwood_box
              LABELS porousmediumflow 1p2c 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=HortonRogersLapwood
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_hortonrogerslapwood_box
              CMD_ARGS params.input -Problem.Name test_1p2c_hortonrogerslapwood_box
                       -Problem.GrowthRateTolerance 1e-2)

# decaying perturbation at Ra = 2 pi^2
dumux_add_test(NAME test_1p2c_hortonrogerslapwood_decay_box
              LABELS porousmediumflow 1p2c 1pnc
              TARGET test_1p2c_hortonrogerslapwood_box
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_hortonrogerslapwood_box
              CMD_ARGS params.input -Problem.Name test_1p2c_hortonrogerslapwood_decay_box
                       -SpatialParams.Permeability 1.6097214109829741e-15
                       -Problem.GrowthRateTolerance 1e-2)

# solution at rest with Scheidegger's dispersion: the velocity and with it the dispersion vanish
dumux_add_test(NAME test_1p2c_hortonrogerslapwood_dispersion_box
              LABELS porousmediumflow 1p2c 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=HortonRogersLapwoodDispersion
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_hortonrogerslapwood_dispersion_box
              CMD_ARGS params.input -Problem.Name test_1p2c_hortonrogerslapwood_dispersion_box
                       -Problem.Amplitude 0
                       -SpatialParams.LongitudinalDispersivity 0.1
                       -SpatialParams.TransverseDispersivity 0.01
                       -Problem.MaxDeviationTolerance 1e-12)
