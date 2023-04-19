# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(1p_1p)
add_subdirectory(1p_richards)
add_subdirectory(1p2c_richards2c)
add_subdirectory(root_soil_benchmark)

dune_symlink_to_source_files(FILES "grids")
