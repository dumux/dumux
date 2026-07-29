# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

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
