# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "sphere.raw")
dumux_add_test(NAME test_ff_stokes_permeability_upscaling_staggered
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=PoreFlowTestMomentumStaggered
               LABELS freeflow navierstokes
               CMAKE_GUARD "( dune-subgrid_FOUND )")

dumux_add_test(NAME test_ff_stokes_permeability_upscaling_diamond
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=PoreFlowTestMomentumDiamond
               LABELS freeflow navierstokes
               CMAKE_GUARD "( dune-subgrid_FOUND )")
