# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input run_benchmark.py linear_stability.py)

# amplification of a finger mode at Ra = 500, k = 8 pi from t = 0.01 to t = 2, reference from linear_stability.py
dumux_add_test(NAME test_1p2c_onesidedconvection_box
              LABELS porousmediumflow 1p2c 1pnc
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_onesidedconvection_box
              CMD_ARGS params.input -Problem.Name test_1p2c_onesidedconvection_box
                       -Problem.ReferenceLogAmplification 1.35860074
                       -Problem.LogAmplificationTolerance 0.025)

# diffusive boundary layer without perturbation: horizontally uniform, flux 1/sqrt(pi Ra t)
dumux_add_test(NAME test_1p2c_onesidedconvection_diffusion_box
              LABELS porousmediumflow 1p2c 1pnc
              TARGET test_1p2c_onesidedconvection_box
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_onesidedconvection_box
              CMD_ARGS params.input -Problem.Name test_1p2c_onesidedconvection_diffusion_box
                       -Problem.Amplitude 0
                       -TimeLoop.DimensionlessEndTime 1
                       -Problem.FluxTolerance 1e-3
                       -Problem.MaxRowDeviationTolerance 1e-10)
