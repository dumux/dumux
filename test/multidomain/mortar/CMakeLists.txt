# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_mortar_trace
              LABELS multidomain multidomain_mortar
              SOURCES test_trace.cc
              CMAKE_GUARD dune-foamgrid_FOUND)

dumux_add_test(NAME test_mortar_decomposition
              LABELS multidomain multidomain_mortar
              SOURCES test_decomposition.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
