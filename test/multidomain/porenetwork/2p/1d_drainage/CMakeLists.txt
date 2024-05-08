# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES extractresults.py porebypore.py eff_acc.py l2error.py)

dumux_add_test(NAME test_pnm_2p_1d_drainage
              SOURCES main.cc
              LABELS porenetwork
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d_drainage
              CMD_ARGS params.input -Pnm.Problem.Name test_pnm_2p_drainage)

dumux_add_test(NAME test_pnm_2p_1d_drainage_reg
              SOURCES main.cc
              LABELS porenetwork
              COMPILE_DEFINITIONS USETHETAREGULARIZATION=1
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d_drainage_reg
              CMD_ARGS params.input -Pnm.Problem.Name test_pnm_2p_drainage_reg)

dumux_add_test(NAME test_pnm_2p_1d_drainage_md
              SOURCES main_md.cc
              LABELS porenetwork
              COMPILE_DEFINITIONS USETHETAREGULARIZATION=1 THROATCONSTRAINT=1
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_1d_drainage_md
              CMD_ARGS params_md.input)
