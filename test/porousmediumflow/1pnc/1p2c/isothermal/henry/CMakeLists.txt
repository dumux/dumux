# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_case2.input)

# The Henry problem benchmark of Fahs et al. (2016, WRR, doi:10.1002/2016WR019288),
# "Test Case 1": the classic (purely diffusive) Henry (1964) configuration. This is
# validated not against our own accepted output but against an independently published
# semianalytical solution (Appendix D, Table D1): validate_fahs2016.py extracts the
# 10/50/90% isochlor positions from our simulation's output and compares them directly
# to the digitized literature values (fahs2016_table_d1.csv). The run is extended well
# past the ~4 h physical time the paper reports for this case to reach steady state, to
# give a safe margin.
dumux_add_test(NAME test_1p2c_henry_fahs_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=HenryFahsTest
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/validate_fahs2016.py
              CMD_ARGS  --vtu ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_box-02000.vtu
                        --reference ${CMAKE_CURRENT_SOURCE_DIR}/fahs2016_table_d1.csv
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_box params.input")

# "Test Case 2": same problem, but with velocity-dependent (Scheidegger) dispersion
# enabled (alphaL=0.1 m, alphaT=0.01 m, smaller molecular diffusion 9.43e-8 m^2/s --
# see params_case2.input and properties.hh's HenryFahsCase2Test). Validated against
# Appendix D, Table D2 the same way Test Case 1 is validated against Table D1. Uses a
# finer mesh and longer run than Test Case 1 (see params_case2.input) since the paper
# reports both mesh sensitivity and time-to-steady-state increase substantially once
# dispersion is enabled.
dumux_add_test(NAME test_1p2c_henry_fahs_case2_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              LABELS 1pnc
              COMPILE_DEFINITIONS TYPETAG=HenryFahsCase2Test
              TIMEOUT 7200 # [s]; finer mesh + 2x the time steps of Test Case 1 (which
                           # took ~117s wall clock), plus dispersion's extra per-flux
                           # cost -- default ctest timeout (1500s) is unlikely to be
                           # enough. Not empirically measured, adjust once you've timed
                           # an actual run.
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/validate_fahs2016.py
              CMD_ARGS  --vtu ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case2_box-04000.vtu
                        --reference ${CMAKE_CURRENT_SOURCE_DIR}/fahs2016_table_d2.csv
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_fahs_case2_box params_case2.input")
