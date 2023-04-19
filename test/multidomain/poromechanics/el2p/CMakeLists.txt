# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_md_poromechanics_el2p
              LABELS multidomain poromechanics 2p poroelastic
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 3600
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el2p_2p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p_twop-00010.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el2p_poroelastic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p_poroelastic-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p params.input
                                                              -Vtk.OutputName test_md_poromechanics_el2p")

dumux_add_test(NAME test_md_poromechanics_el2p_gravity
              LABELS multidomain poromechanics 2p poroelastic
              TARGET test_md_poromechanics_el2p
              CMAKE_GUARD HAVE_UMFPACK
              TIMEOUT 3600
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el2p_gravity_2p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p_gravity_twop-00010.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_poromechanics_el2p_gravity_poroelastic-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p_gravity_poroelastic-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_poromechanics_el2p params.input
                                                              -Vtk.OutputName test_md_poromechanics_el2p_gravity -Problem.EnableGravity true -Newton.MaxRelativeShift 1e-8")
