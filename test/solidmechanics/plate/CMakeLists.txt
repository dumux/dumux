# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES convergencetest.py plot_convergence.py regenerate_doc_images.py disk.geo)

add_subdirectory(kirchhoff_love)
add_subdirectory(membrane)
add_subdirectory(mindlin_reissner)
