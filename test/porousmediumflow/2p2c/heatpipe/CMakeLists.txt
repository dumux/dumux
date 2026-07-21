# SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_heatpipe_box
              TIMEOUT 1800
              SOURCES heatpipe.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_heatpipe_box
              CMD_ARGS heatpipe.input -TimeLoop.TEnd 1e5 -TimeLoop.EpisodeLength 1e5)
