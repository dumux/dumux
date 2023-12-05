# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_md_boundary_ff1pni_pnm1pni
               SOURCES main.cc
               COMPILE_DEFINITIONS ISOTHERMAL=0) #TODO: add reference files as soon as solutions seem reasonable
