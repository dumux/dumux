# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(angeli)
add_subdirectory(channel)
add_subdirectory(donea)
add_subdirectory(kovasznay)
add_subdirectory(periodic)
add_subdirectory(sincos)
add_subdirectory(unstructured)

dune_symlink_to_source_files(FILES convergence.sh)
