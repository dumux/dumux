# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_3d.input)

dumux_add_test(NAME test_hyperelastic_growth_layered_shell_box_2d
               LABELS solidmechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=LayeredShellBox DIMENSION=2
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../convergencetest.py
               CMD_ARGS test_hyperelastic_growth_layered_shell_box_2d params.input "6 4 16" 1.9 0.9)

dumux_add_test(NAME test_hyperelastic_growth_layered_shell_pq2hybrid_2d
               LABELS solidmechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=LayeredShellPQ2Hybrid DIMENSION=2
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../convergencetest.py
               CMD_ARGS test_hyperelastic_growth_layered_shell_pq2hybrid_2d params.input "6 4 16" 1.8 1.4)

dumux_add_test(NAME test_hyperelastic_growth_layered_shell_box_3d
               LABELS solidmechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=LayeredShellBox DIMENSION=3
               CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../convergencetest.py
               CMD_ARGS test_hyperelastic_growth_layered_shell_box_3d params_3d.input "3 2 4" 1.8 0.9)
