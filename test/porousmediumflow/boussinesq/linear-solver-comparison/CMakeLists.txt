# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

# Not a physics/regression test: benchmarks several ISTL linear solvers
# against matrices assembled during a transient Boussinesq run (both
# formulations). Split into two stages:
#
#   1. "generate" runs the transient simulation (see white-noise-perturbations
#      for the same physical setup) to TimeLoop.TEnd, exporting the assembled
#      Jacobian/residual at each time in MatrixExport.Times to MatrixMarket
#      files (matrixio.hh).
#   2. "solve" reads those files back in and benchmarks each candidate solver
#      against them, reporting solver iterations/reduction/wall time to
#      <name>_linsolver_comparison.csv (solvercompare.hh).
#
# The solve tests DEPENDS on their matching generate test having run first
# so the matrix files exist (same ordering mechanism as
# test_boussinesq_dissolution_comparison_eval in the parent CMakeLists.txt).

# --- pressure formulation -----------------------------------------------

dumux_add_test(NAME test_boussinesq_linsolver_generate_pressure
              LABELS porousmediumflow boussinesq
              SOURCES main_generate_pressure.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_linsolver_generate_pressure
              CMD_ARGS params.input -Problem.Name test_boussinesq_linsolver_compare_pressure)

dumux_add_test(NAME test_boussinesq_linsolver_solve_pressure
              LABELS porousmediumflow boussinesq
              SOURCES main_solve_pressure.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_linsolver_solve_pressure
              CMD_ARGS params.input -Problem.Name test_boussinesq_linsolver_compare_pressure)

set_tests_properties(test_boussinesq_linsolver_solve_pressure
    PROPERTIES DEPENDS test_boussinesq_linsolver_generate_pressure)

# --- vorticity formulation -----------------------------------------------

dumux_add_test(NAME test_boussinesq_linsolver_generate_vorticity
              LABELS porousmediumflow boussinesq
              SOURCES main_generate_vorticity.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_linsolver_generate_vorticity
              CMD_ARGS params.input -Problem.Name test_boussinesq_linsolver_compare_vorticity)

dumux_add_test(NAME test_boussinesq_linsolver_solve_vorticity
              LABELS porousmediumflow boussinesq
              SOURCES main_solve_vorticity.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_linsolver_solve_vorticity
              CMD_ARGS params.input -Problem.Name test_boussinesq_linsolver_compare_vorticity)

set_tests_properties(test_boussinesq_linsolver_solve_vorticity
    PROPERTIES DEPENDS test_boussinesq_linsolver_generate_vorticity)
