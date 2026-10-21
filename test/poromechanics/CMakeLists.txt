# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES sphere_eighth.geo sphere_eighth_coarse.msh)

add_subdirectory(poroelastic)
add_subdirectory(poroelastic_large_deformations)
