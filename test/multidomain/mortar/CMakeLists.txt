# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_mortar_decomposition
              LABELS multidomain multidomain_mortar
              SOURCES test_decomposition.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
dumux_add_test(NAME test_mortar_model
              LABELS multidomain multidomain_mortar
              SOURCES test_model.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
