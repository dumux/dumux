# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "X-dir-CapillaryPressureVsAverageSaturation.dat")
dune_symlink_to_source_files(FILES "X-dir-Non-wetting-EffPermeabilityVsAverageSaturation.dat")
dune_symlink_to_source_files(FILES "X-dir-Wetting-EffPermeabilityVsAverageSaturation.dat")

# using box
dumux_add_test(NAME test_2p_incompressible_box_printmed
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleBox
              CMAKE_GUARD dune-subgrid_FOUND
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box params.input -Problem.Name test_2p_incompressible_box")

dumux_add_test(NAME test_2p_incompressible_box_analytic_printmed
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPIncompressibleBox DIFFMETHOD=DiffMethod::analytic
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box_analytic-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box_analytic params.input -Problem.Name test_2p_incompressible_box_analytic -Newton.EnablePartialReassembly false")
