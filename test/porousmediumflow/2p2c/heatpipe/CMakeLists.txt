# SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_heatpipe_box
              TIMEOUT 1800
              SOURCES heatpipe.cc
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/lecture/references/heatpipe-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/heatpipe-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/heatpipe -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/heatpipe.input -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/heatpipe.dgf -TimeLoop.TEnd 1e5 -TimeLoop.EpisodeLength  1e5")
