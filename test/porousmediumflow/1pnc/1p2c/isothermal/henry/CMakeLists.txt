# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_case2.input)

# The Henry problem benchmark of Fahs et al. (2016, WRR, doi:10.1002/2016WR019288),
# "Test Case 1": the classic (purely diffusive) Henry (1964) configuration. Built here
# (SOURCES) but *not* run directly: the actual simulation run lives on the regression
# target below, which this target depends on and reuses the VTU output of (see DEPENDS
# below), avoiding a second full simulation run of the same case. This target instead
# validates that already-produced output against an independently published
# semianalytical solution (Appendix D, Table D1): validate_fahs2016.py extracts the
# 10/50/90% isochlor positions from the simulation output and compares them directly
# to the digitized literature values (fahs2016_table_d1.csv). The run is extended well
# past the ~4 h physical time the paper reports for this case to reach steady state, to
# give a safe margin.
dumux_add_test(NAME test_1p2c_henry_fahs_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=HenryFahsTest
              TIMEOUT 60 # [s]; no simulation run here anymore, just reading an
                         # already-produced VTU and interpolating a handful of isochlors
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/validate_fahs2016.py
              CMD_ARGS  --vtu ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_box-00100.vtu
                        --reference ${CMAKE_CURRENT_SOURCE_DIR}/fahs2016_table_d1.csv)

# Regression test for Test Case 1: reuses the executable built above (no SOURCES, so
# no rebuild), runs it once under ${DUMUX_RUNTEST}'s standard `--script fuzzy`
# mechanism (see e.g. test/porousmediumflow/co2/CMakeLists.txt, test/porousmediumflow/
# 2pncmin/isothermal/CMakeLists.txt for the same pattern), comparing the full VTU
# output field-by-field against a stored, accepted reference
# (test/references/test_1p2c_henry_fahs_box-reference.vtu, generated at the current
# 240x80/1 d setup). This is the one and only run of the simulation for this test case
# -- test_1p2c_henry_fahs_box above (DEPENDS on this target) reuses its VTU output
# rather than re-running, since validate_fahs2016.py only checks a handful of
# extracted isochlor positions, not the full solution field that this target checks.
dumux_add_test(NAME test_1p2c_henry_fahs_box_regression
              TARGET  test_1p2c_henry_fahs_box
              LABELS porousmediumflow 1p2c 1pnc
              TIMEOUT 7200 # [s]; measured ~953s wall clock with the direct (UMFPack)
                           # linear solver on the current 240x80/1 d setup
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_henry_fahs_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_box-00100.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_box params.input")
set_tests_properties(test_1p2c_henry_fahs_box PROPERTIES DEPENDS test_1p2c_henry_fahs_box_regression)

# "Test Case 2": same problem, but with velocity-dependent (Scheidegger) dispersion
# enabled (alphaL=0.1 m, alphaT=0.01 m, smaller molecular diffusion 9.43e-8 m^2/s --
# see params_case2.input and properties.hh's HenryFahsCase2Test). Validated against
# Appendix D, Table D2 the same way Test Case 1 is validated against Table D1. Same
# mesh and number of time steps as Test Case 1 (see params_case2.input). Built here
# (SOURCES) but *not* run directly, same DEPENDS-on-the-regression-target pattern as
# test_1p2c_henry_fahs_box above -- see that target's comment.
dumux_add_test(NAME test_1p2c_henry_fahs_case2_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=HenryFahsCase2Test
              TIMEOUT 60 # [s]; no simulation run here anymore, see test_1p2c_henry_fahs_box
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/validate_fahs2016.py
              CMD_ARGS  --vtu ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case2_box-00100.vtu
                        --reference ${CMAKE_CURRENT_SOURCE_DIR}/fahs2016_table_d2.csv)

# Regression test for Test Case 2, same pattern as test_1p2c_henry_fahs_box_regression
# above, against test/references/test_1p2c_henry_fahs_case2_box-reference.vtu. This is
# the one and only run of the simulation for this test case; dispersion's extra
# per-flux cost makes this run somewhat slower than Test Case 1's (measured: ~1269s vs.
# Test Case 1's ~953s).
dumux_add_test(NAME test_1p2c_henry_fahs_case2_box_regression
              TARGET  test_1p2c_henry_fahs_case2_box
              LABELS porousmediumflow 1p2c 1pnc
              TIMEOUT 7200 # [s]; measured ~1269s wall clock with the direct (UMFPack)
                           # linear solver on the current 240x80/1 d setup
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_henry_fahs_case2_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case2_box-00100.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case2_box params_case2.input")
set_tests_properties(test_1p2c_henry_fahs_case2_box PROPERTIES DEPENDS test_1p2c_henry_fahs_case2_box_regression)