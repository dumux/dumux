# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(boundaryfluxes)
add_subdirectory(bowl)
add_subdirectory(dambreak)
add_subdirectory(fluxlimiter)
add_subdirectory(iofields)
add_subdirectory(jacobian)
add_subdirectory(lakeatrest)
add_subdirectory(manufactured)
add_subdirectory(poiseuilleflow)
add_subdirectory(reconstruction)
add_subdirectory(roughchannel)

dune_symlink_to_source_files(FILES plotstyle.py)
