# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

dumux_add_test(
    NAME test_gripper
    LABELS solidmechanics poromechanics
    SOURCES main.cc
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_gripper
    CMD_ARGS params.input
)

#target_compile_options(test_gripper PRIVATE -fsanitize=address)
#target_link_options(test_gripper PRIVATE -fsanitize=address)
