# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME utrecht_exp_14_fin
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pnc-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc-00010.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_exp_14_fir
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 USETHETAREGULARIZATION=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pnc-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc-00010.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_exp_14_md
               SOURCES main_md.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1 USETHETAREGULARIZATION=1 THROATCONSTRAINT=1
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pnc-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc-00010.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc")

dumux_add_test(NAME utrecht_exp_14_ni
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2pnc_ni-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc_ni-00004.vtp
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2pnc_ni params_ni.input -Problem.Name test_pnm_2pnc_ni")
