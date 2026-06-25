# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# streamfunction / vector-potential vorticity formulation
add_subdirectory(vorticity)

# pressure-based formulation with Nitsche pressure fix
add_subdirectory(pressure)

# ----- Comparison tests (Ra=8000, 200×200) -----------------------------------
# Step 1a: run vorticity formulation
dumux_add_test(NAME test_boussinesq_comparison_run_vorticity
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_onesided_rb_box
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/vorticity/test_boussinesq_onesided_rb_box
              CMD_ARGS params.input -Problem.Name comparison_vorticity)

set_tests_properties(test_boussinesq_comparison_run_vorticity
    PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vorticity)

# Step 1b: run pressure formulation
dumux_add_test(NAME test_boussinesq_comparison_run_pressure
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_pressure_box
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/pressure/test_boussinesq_pressure_box
              CMD_ARGS params.input -Problem.Name comparison_pressure)

set_tests_properties(test_boussinesq_comparison_run_pressure
    PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pressure)

# Step 2: read the output and produce comparison figures
# Requires both step 1a and 1b to have completed
dumux_add_test(NAME test_boussinesq_comparison_eval
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_onesided_rb_box
              COMMAND ${Python3_EXECUTABLE}
              CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/eval_comparison.py
                       ${CMAKE_CURRENT_BINARY_DIR}
                       --out ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(test_boussinesq_comparison_eval
    PROPERTIES DEPENDS "test_boussinesq_comparison_run_vorticity;test_boussinesq_comparison_run_pressure")
