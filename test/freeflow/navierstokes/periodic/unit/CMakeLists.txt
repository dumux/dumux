# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "periodic.dgf")
dumux_add_test(NAME test_ff_navierstokes_periodic
               LABELS freeflow navierstokes
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-spgrid_FOUND )"
               COMPILE_DEFINITIONS UPWINDSCHEMEORDER=1
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --zeroThreshold {"velocity_liq \(m/s\)":1e-12}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_periodic-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_periodic-00000.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_periodic params.input
                                        -Problem.Name test_ff_navierstokes_periodic")

add_executable(test_ff_navierstokes_periodic_subgrid EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_navierstokes_periodic_subgrid PUBLIC UPWINDSCHEMEORDER=1 USESUBGRID=1)

dumux_add_test(NAME test_ff_navierstokes_periodic_subgrid_trivial
               TARGET test_ff_navierstokes_periodic_subgrid
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-spgrid_FOUND AND dune-subgrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --zeroThreshold {"velocity_liq \(m/s\)":1e-12}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_periodic-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_periodic_subgrid_trivial-00000.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_periodic_subgrid params.input
                                        -Problem.Name test_ff_navierstokes_periodic_subgrid_trivial
                                        -Problem.TrivialSelector true")

dumux_add_test(NAME test_ff_navierstokes_periodic_subgrid_nontrivial
               TARGET test_ff_navierstokes_periodic_subgrid
               LABELS freeflow navierstokes
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-spgrid_FOUND AND dune-subgrid_FOUND )"
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS      --script fuzzy --zeroThreshold {"velocity_liq \(m/s\)":1e-12}
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_periodic_subgrid-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_periodic_subgrid-00000.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_periodic_subgrid params.input
                                        -Problem.Name test_ff_navierstokes_periodic_subgrid
                                        -Problem.TrivialSelector false")
