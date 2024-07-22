# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES generate_lattice_grid.py plot_sn_t.py 1d_4pores.dgf)

dumux_add_test(NAME pnm_3d_drainage
               SOURCES main.cc
               LABELS porenetwork
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/pnm_3d params.input")

dumux_add_test(NAME pnm_3d_drainage_reg
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS USETHETAREGULARIZATION=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/pnm_3d_reg params.input")

dumux_add_test(NAME pnm_3d_drainage_md
               SOURCES main_md.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS USETHETAREGULARIZATION=1 THROATCONSTRAINT=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/pnm_3d_md params_md.input")
