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

# ----- Dissolution-flux comparison tests (Ra=8000, parallel + h-adaptive + rebalanced) ---
# Both formulations now run on ALUGrid with the same parallel/h-adaptive/rebalancing setup
# as their production test targets (see white-noise-perturbations/CMakeLists.txt) -- fixed
# at MPI_RANKS 4 (a single value, not a list: dune_add_test appends a "-mpi-<n>" suffix to
# the actual ctest name per rank count, and every rank count would otherwise share the same
# -Problem.Name and so the same output files, silently clobbering each other if more than
# one were run). Needs DUNE_MAX_TEST_CORES >= 4 (dune-common default is 2) for this rank
# count to actually register as a CTest entry -- see cmake.opts / your CMake profile.
# A serial run never has to reconcile ILU's blockwise-per-rank conditioning or an actual
# grid repartition, so it doesn't exercise the same code paths as production use -- these
# metrics are computed under the real configuration instead.
#
# Step 1a: run vorticity formulation
dumux_add_test(NAME test_boussinesq_dissolution_comparison_run_vorticity
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_vorticity_box
              CMAKE_GUARD dune-alugrid_FOUND
              MPI_RANKS 4
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations/test_boussinesq_white_noise_vorticity_box
              CMD_ARGS params.input -Problem.Name comparison_vorticity
                       -Adaptive.EnableDynamicLoadBalancing true
                       -Adaptive.RebalanceEvery 3
                       -Adaptive.VerticalStripPartitioning true)

set_tests_properties(test_boussinesq_dissolution_comparison_run_vorticity-mpi-4
    PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations)

# Step 1b: run pressure formulation
# Uses the _amg executable, not the default ILU one: ILU+GMRES(30) hits a hard convergence
# wall on this operator at Ra=8000 (confirmed -- Newton never recovers, halves dt to the
# noise floor and aborts); AMG+GMRES (test_boussinesq_white_noise_pressure_box_amg, see
# white-noise-perturbations/CMakeLists.txt) converges cleanly. Same fix as proven for the
# cosine-perturbations pressure test at Ra=500.
dumux_add_test(NAME test_boussinesq_dissolution_comparison_run_pressure
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_pressure_box_amg
              CMAKE_GUARD dune-alugrid_FOUND
              MPI_RANKS 4
              TIMEOUT 200000
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations/test_boussinesq_white_noise_pressure_box_amg
              CMD_ARGS params.input -Problem.Name comparison_pressure
                       -Adaptive.EnableDynamicLoadBalancing true
                       -Adaptive.RebalanceEvery 3
                       -Adaptive.VerticalStripPartitioning true)

set_tests_properties(test_boussinesq_dissolution_comparison_run_pressure-mpi-4
    PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/white-noise-perturbations)

# Step 2: read the output and produce comparison figures
# Requires both step 1a and 1b to have completed
dumux_add_test(NAME test_boussinesq_dissolution_comparison_eval
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_white_noise_vorticity_box
              CMAKE_GUARD dune-alugrid_FOUND
              COMMAND ${Python3_EXECUTABLE}
              CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/eval_comparison.py
                       ${CMAKE_CURRENT_BINARY_DIR}
                       --out ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(test_boussinesq_dissolution_comparison_eval
    PROPERTIES DEPENDS "test_boussinesq_dissolution_comparison_run_vorticity-mpi-4;test_boussinesq_dissolution_comparison_run_pressure-mpi-4")
