# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_richards_annulus
               LABELS porousmediumflow richards
               SOURCES main.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_richards_annulus
               CMD_ARGS params.input -Problem.Name test_richards_annulus -SpatialParams.SoilType Clay)

dumux_add_test(NAME test_richards_annulus_loam
               LABELS porousmediumflow richards
               TARGET test_richards_annulus
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_richards_annulus
               CMD_ARGS params.input -Problem.Name test_richards_annulus -SpatialParams.SoilType Loam)

dumux_add_test(NAME test_richards_annulus_instationary
               LABELS porousmediumflow richards
               SOURCES main_instationary.cc
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_richards_annulus_instationary
               CMD_ARGS params.input -Problem.Name test_richards_annulus_instationary -SpatialParams.SoilType Clay)

dumux_add_test(NAME test_richards_annulus_instationary_loam
               LABELS porousmediumflow richards
               TARGET test_richards_annulus_instationary
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_richards_annulus_instationary
               CMD_ARGS params.input -Problem.Name test_richards_annulus_instationary -SpatialParams.SoilType Loam)
