# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# pressure/, vorticity/ and common/ hold only shared, model-related headers
# (Darcy law, problem.hh, model code, fluidsystem, spatialparams) -- no
# CMakeLists.txt / build targets, no main.cc, of their own. Every executable
# (main.cc) lives in one of the application dirs below, one per formulation.

# dissolution-flux comparison (amplitude=0 reduces to the unperturbed erfc
# baseline) and Riaz et al. (2006) finger-count/tip validation, both formulations
add_subdirectory(white-noise-perturbations)

# single-sinusoidal-mode-perturbed growth-rate validation (Riaz et al. 2006
# Fig. 8), both formulations
add_subdirectory(cosine-perturbations)

# linear-solver benchmark (small grid, both formulations) -- not a physics test
add_subdirectory(linear-solver-comparison)

# ----- Dissolution-flux comparison tests (Ra=8000, 200×200) ------------------
# Step 1a: run vorticity formulation
dumux_add_test(NAME test_boussinesq_dissolution_comparison_run_vorticity
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_vorticity_box
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations/test_boussinesq_white_noise_vorticity_box
              CMD_ARGS params.input -Problem.Name comparison_vorticity)

set_tests_properties(test_boussinesq_dissolution_comparison_run_vorticity
    PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations)

# Step 1b: run pressure formulation
dumux_add_test(NAME test_boussinesq_dissolution_comparison_run_pressure
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_pressure_box
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations/test_boussinesq_white_noise_pressure_box
              CMD_ARGS params.input -Problem.Name comparison_pressure)

set_tests_properties(test_boussinesq_dissolution_comparison_run_pressure
    PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations)

# Step 2: read the output and produce comparison figures
# Requires both step 1a and 1b to have completed
dumux_add_test(NAME test_boussinesq_dissolution_comparison_eval
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_vorticity_box
              COMMAND ${Python3_EXECUTABLE}
              CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/eval_comparison.py
                       ${CMAKE_CURRENT_BINARY_DIR}
                       --out ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(test_boussinesq_dissolution_comparison_eval
    PROPERTIES DEPENDS "test_boussinesq_dissolution_comparison_run_vorticity;test_boussinesq_dissolution_comparison_run_pressure")
