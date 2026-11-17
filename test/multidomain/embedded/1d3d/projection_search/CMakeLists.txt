# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_md_embedded_1d3d_projection_search
               SOURCES main.cc
               LABELS unit multidomain multidomain_embedded
               CMAKE_GUARD dune-foamgrid_FOUND)
