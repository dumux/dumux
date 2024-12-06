# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(darcy_darcy)

dumux_add_test(NAME test_mortar_decomposition
              LABELS multidomain multidomain_mortar
              SOURCES test_decomposition.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
dumux_add_test(NAME test_mortar_model
              LABELS multidomain multidomain_mortar
              SOURCES test_model.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
dumux_add_test(NAME test_mortar_interface_operator
              LABELS multidomain multidomain_mortar
              SOURCES test_interface_operator.cc
              CMAKE_GUARD dune-foamgrid_FOUND)
