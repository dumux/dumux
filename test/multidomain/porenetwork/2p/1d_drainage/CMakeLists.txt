# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES extractresults.py porebypore.py eff_acc_largeDt.py compare_convergence.py)

dumux_add_test(NAME test_pnm_2p_1d_drainage
              SOURCES main.cc
              LABELS multidomain porenetwork
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files  ${CMAKE_SOURCE_DIR}/test/references/test_multidomain_porenetwork_2p_1d_drainage_fin-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_fi-00050.vtp
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d_drainage params.input -Pnm.Problem.Name test_pnm_2p_fi")

dumux_add_test(NAME test_pnm_2p_1d_drainage_reg
              SOURCES main.cc
              LABELS porenetwork
              COMPILE_DEFINITIONS USETHETAREGULARIZATION=1
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files  ${CMAKE_SOURCE_DIR}/test/references/test_multidomain_porenetwork_2p_1d_drainage_fir-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_reg-00344.vtp
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d_drainage_reg params.input -Pnm.Problem.Name test_pnm_2p_reg")

dumux_add_test(NAME test_pnm_2p_1d_drainage_md
              SOURCES main_md.cc
              LABELS porenetwork
              COMPILE_DEFINITIONS USETHETAREGULARIZATION=1 THROATCONSTRAINT=1
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files  ${CMAKE_SOURCE_DIR}/test/references/test_multidomain_porenetwork_2p_1d_drainage_md-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_md-00053.vtp
                                     ${CMAKE_SOURCE_DIR}/test/references/test_multidomain_porenetwork_2p_1d_drainage_md_constraint-reference.vtp
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_md_constraint-00053.vtp
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d_drainage_md params_md.input")
