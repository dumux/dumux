# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
dumux_add_test(NAME test_ff_navierstokes_taylorcouette
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG_MOMENTUM=TaylorCouetteTestMomentumDiamond
               COMPILE_DEFINITIONS TYPETAG_MASS=TaylorCouetteTestMassTpfa
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-uggrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_taylorcouette.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_taylorcouette-00001.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_taylorcouette params.input
                             -Problem.Name test_ff_navierstokes_taylorcouette -Problem.EnableInertiaTerms true")
dune_symlink_to_source_files(FILES "params.input")