# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_md_boundary_ff1p_pnm1p_3D
               SOURCES main.cc
               LABELS multidomain multidomain_boundary freeflowpnm 1p navierstokes porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
        )
