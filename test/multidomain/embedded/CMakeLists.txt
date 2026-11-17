# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(1d3d)
add_subdirectory(2d3d)

dumux_add_test(NAME test_integrationpointsource
              SOURCES test_integrationpointsource.cc
              LABELS unit multidomain multidomain_embedded)
