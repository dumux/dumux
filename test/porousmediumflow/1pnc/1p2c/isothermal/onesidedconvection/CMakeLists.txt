# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input linear_stability.py)

# growth of a perturbation of the boundary layer at Ra = 500, k = 8 pi from t = 0.01 to t = 2, reference from linear_stability.py
dumux_add_test(NAME test_1p2c_onesidedconvection_box
              LABELS porousmediumflow 1p2c 1pnc
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_onesidedconvection_box
              CMD_ARGS params.input -Problem.Name test_1p2c_onesidedconvection_box
                       -Problem.ReferenceLogAmplification 1.35860074
                       -Problem.LogAmplificationTolerance 0.025)
