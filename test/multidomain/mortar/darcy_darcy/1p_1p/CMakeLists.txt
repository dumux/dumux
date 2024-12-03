# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_input_file_links()
dumux_add_test(NAME test_mortar_darcy_1p
              LABELS multidomain multidomain_mortar
              SOURCES main.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
