# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
dune_symlink_to_source_files(FILES params.input visualize_wood.py)
dumux_add_test(
    NAME test_wood
    LABELS solidmechanics elastic
    SOURCES main.cc
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_wood
    CMD_ARGS params.input -TimeLoop.TEnd 1e6
)
