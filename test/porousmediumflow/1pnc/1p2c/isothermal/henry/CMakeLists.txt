# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params_case1.input)

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
dumux_add_test(NAME test_1p2c_henry_fahs_case1_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=HenryFahsTest
              TIMEOUT 60 # [s]; no simulation run here anymore, just reading an
                         # already-produced VTU and interpolating a handful of isochlors
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/validate_fahs2016.py
              CMD_ARGS  --vtu ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case1_box-00100.vtu
                        --reference ${CMAKE_CURRENT_SOURCE_DIR}/fahs2016_table_d1.csv)

# Regression test for Test Case 1: reuses the executable built above (no SOURCES, so
# no rebuild), runs it once under ${DUMUX_RUNTEST}'s standard `--script fuzzy`
# mechanism (see e.g. test/porousmediumflow/co2/CMakeLists.txt, test/porousmediumflow/
# 2pncmin/isothermal/CMakeLists.txt for the same pattern), comparing the full VTU
# output field-by-field against a stored, accepted reference
# (test/references/test_1p2c_henry_fahs_case1_box-reference.vtu, generated at the current
# 120x40/1 d setup). This is the one and only run of the simulation for this test case
# -- test_1p2c_henry_fahs_case1_box above (DEPENDS on this target) reuses its VTU output
# rather than re-running, since validate_fahs2016.py only checks a handful of
# extracted isochlor positions, not the full solution field that this target checks.
dumux_add_test(NAME test_1p2c_henry_fahs_case1_box_regression
              TARGET  test_1p2c_henry_fahs_case1_box
              LABELS porousmediumflow 1p2c 1pnc
              TIMEOUT 7200 # [s]; measured ~953s wall clock with the direct (UMFPack)
                           # linear solver on the earlier 240x80/1 d setup
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2c_henry_fahs_case1_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case1_box-00100.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case1_box params_case1.input")
set_tests_properties(test_1p2c_henry_fahs_case1_box PROPERTIES DEPENDS test_1p2c_henry_fahs_case1_box_regression)
