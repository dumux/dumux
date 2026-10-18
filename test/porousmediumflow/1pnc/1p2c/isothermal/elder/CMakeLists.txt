# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input run_benchmark.py elder_pseudospectral.py)

# Ra = 60: solute mass and Sherwood number at t = 0.03, 0.1, 0.3 against the pseudo-spectral reference solution
dumux_add_test(NAME test_1p2c_elder_box
              LABELS porousmediumflow 1p2c 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=Elder
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_elder_box
              CMD_ARGS params.input -Problem.Name test_1p2c_elder_box
                       -Problem.ReferenceMass "0.2056540334 0.4180151245 0.6138850339"
                       -Problem.ReferenceSherwood "1.911360151 1.529963756 1.252857758"
                       -Problem.RelativeTolerance 5e-3)

# Ra = 0: relative L2 error of the mass fraction against the exact solution
dumux_add_test(NAME test_1p2c_elder_nobuoyancy_box
              LABELS porousmediumflow 1p2c 1pnc
              TARGET test_1p2c_elder_box
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_elder_box
              CMD_ARGS params.input -Problem.Name test_1p2c_elder_nobuoyancy_box
                       -Problem.DensityDifference 0
                       -Problem.ExactSolutionTolerance 1.5e-3)
