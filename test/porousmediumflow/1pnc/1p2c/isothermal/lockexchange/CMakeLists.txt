# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

# full (variable-density) model: density = rho0*(1+C) everywhere, standard CVFE Darcy
# law. Regular regression test against its own accepted reference solution.
dumux_add_test(NAME test_1p2c_lockexchange_full_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=LockExchangeFullTest
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_lockexchange_full_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_full_box-00100.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_full_box params.input -Problem.Name test_1p2c_lockexchange_full_box")

# Boussinesq-approximated variant of the same problem (same mesh/BCs/IC/time steps):
# AdvectionType = BoussinesqCVFEDarcyLaw with a constant reference density, buoyancy
# driven by Problem::buoyantDensity(). Regular regression test against its own
# accepted reference solution, same as the full model above -- this only catches
# unintended changes to the Boussinesq code path itself, it says nothing about how
# well the approximation matches the full model (see the comparison test below for
# that).
dumux_add_test(NAME test_1p2c_lockexchange_boussinesq_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=LockExchangeBoussinesqTest
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_lockexchange_boussinesq_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_boussinesq_box-00100.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_boussinesq_box params.input -Problem.Name test_1p2c_lockexchange_boussinesq_box")

# run after the two reference-based tests above, so all three appear in that order
# under `ctest -R lockexchange`
set_tests_properties(test_1p2c_lockexchange_boussinesq_box
                      PROPERTIES DEPENDS test_1p2c_lockexchange_full_box)

# -----------------------------------------------------------------------------------
# Boussinesq vs. full model: NOT a regression test in the usual sense. The two are
# different models by construction (one keeps density constant outside the buoyancy
# term, the other doesn't), so an exact/tight match is neither expected nor the point.
# This check instead bounds how far the Boussinesq approximation is allowed to drift
# from the full model, at this test's deliberately small (1%, see problem.hh) density
# variance -- i.e. it tests the *quality of the approximation*, not correctness of
# either model in isolation. A future change that grows this gap (e.g. by increasing
# the density variance, or introducing an actual bug in buoyantDensity()) should show
# up here even if both models still pass their own reference checks above.
#
# "rho" is excluded from the comparison: it is *supposed* to differ between the two
# models -- in the Boussinesq variant the tracer-driven part of the density has been
# moved out of the fluid system (volVars.density(), which is what the "rho" VTK field
# is) and into Problem::buoyantDensity(), which only re-enters through the flux, not
# through volVars.density() itself.
#
# --relative is deliberately looser (1.5e-2) than the default (1e-2) used by the
# reference-based tests above: even at this test's small (1%, see problem.hh) density
# variance, the full and Boussinesq models differ somewhat (empirically, up to ~1.3%
# in the tracer fraction fields, at ordinary, non-near-zero magnitudes), which is
# expected front-position sensitivity between two slightly different formulations, not
# a bug. 1.5e-2 is the tightest value (bisected empirically against this scenario's
# actual output) that still comfortably clears that ~1.3% with a small margin for
# platform/compiler floating-point differences -- tightened down from an initial,
# untuned 3e-2 once the approximation was confirmed to behave as expected.
#
# zeroThreshold gives specific fields an absolute floor below which deviations are not
# counted as a failure at all -- i.e. for values at or below the floor, the check is
# effectively absolute rather than relative. This matters for the velocity field's
# minor (cross-flow) component: it stays small throughout (order 1e-7 m/s, versus
# ~1e-6 m/s for the dominant along-flow component), so a purely relative check blows
# up small absolute differences there into large percentages. Both floors below are
# likewise bisected to the tightest value that still passes (2e-8 for velocity, 1e-5
# for the tracer fractions -- the latter is load-bearing, not just conservative:
# tightening it further reintroduces failures from genuinely near-zero trace fraction
# values ahead of the current's leading edge), plus a small margin.
dumux_add_test(NAME test_1p2c_lockexchange_boussinesq_vs_full
              TARGET test_1p2c_lockexchange_boussinesq_box
              LABELS porousmediumflow 1p2c
              LABELS 1pnc
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --relative 1.5e-2
                        --ignore rho
                        --zeroThreshold {"x^solute_liq":1e-5,"X^solute_liq":1e-5,"velocity_liq \(m/s\)":3e-8}
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_full_box-00100.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_boussinesq_box-00100.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_lockexchange_boussinesq_box params.input -Problem.Name test_1p2c_lockexchange_boussinesq_box")

# needs the full model's reference-checked output to already exist; also keeps the
# three tests running in a fixed order (full, then boussinesq, then the comparison)
set_tests_properties(test_1p2c_lockexchange_boussinesq_vs_full
                      PROPERTIES DEPENDS test_1p2c_lockexchange_boussinesq_box)
