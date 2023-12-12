# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES extractresults.py analyze.py)

dumux_add_test(NAME test_pnm_2p_reg
              SOURCES main.cc
              LABELS porenetwork
              COMPILE_DEFINITIONS ISOTHERMAL=1
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p
              CMD_ARGS params.input)
