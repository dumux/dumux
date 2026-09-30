# SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_heatpipe_box
              TIMEOUT 1800
              SOURCES main.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_heatpipe_box
              CMD_ARGS params.input)

# semi-analytical reference solution used by compile_run_plot.py
dumux_add_test(NAME test_heatpipe_odesolver
              SOURCES test_heatpipe_odesolver.cc
              LABELS porousmediumflow 2p2c
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_heatpipe_odesolver
              CMD_ARGS params.input)
