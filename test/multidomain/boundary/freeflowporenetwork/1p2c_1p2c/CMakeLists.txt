# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory("non-isothermal")
add_input_file_links()

dumux_add_test(NAME test_md_boundary_ff1p2c_pnm1p2c
                SOURCES main.cc
                COMPILE_DEFINITIONS ISOTHERMAL=1) #TODO: add reference files as soon as solutions seem reasonable)

dumux_add_test(NAME test_md_boundary_ff1p2c_pnm1p2cni
                SOURCES main.cc
                COMPILE_DEFINITIONS ISOTHERMAL=0) #TODO: add reference files as soon as solutions seem reasonable)
