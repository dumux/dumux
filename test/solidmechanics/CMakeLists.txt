# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES cooksmembrane.msh)

add_subdirectory(elastic)
add_subdirectory(hyperelastic)
add_subdirectory(hyperelastic_voliso)
add_subdirectory(turek_hron)
add_subdirectory(plate)
