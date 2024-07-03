# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_decomposition.cc
               LABELS unit mortar
               CMAKE_GUARD "dune-foamgrid_FOUND")

dumux_add_test(SOURCES test_projectors.cc
               LABELS unit)
