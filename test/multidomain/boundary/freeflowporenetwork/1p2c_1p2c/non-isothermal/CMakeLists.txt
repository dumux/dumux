# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
add_input_file_links()

dumux_add_test(NAME test_md_boundary_ff1p2cni_pnm1p2cni
                SOURCES main.cc) #TODO: add reference files as soon as solutions seem reasonable)
