# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params_diffusion.input params_dispersion.input
                                   fahs2016_table_d1.csv fahs2016_table_d2.csv plot_profiles.py)

dumux_add_test(NAME test_1p2c_henry_diffusion_box
              LABELS porousmediumflow 1p2c 1pnc
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_diffusion_box
              CMD_ARGS params_diffusion.input)

dumux_add_test(NAME test_1p2c_henry_dispersion_box
              LABELS porousmediumflow 1p2c 1pnc
              TARGET test_1p2c_henry_diffusion_box
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_1p2c_henry_diffusion_box
              CMD_ARGS params_dispersion.input)
