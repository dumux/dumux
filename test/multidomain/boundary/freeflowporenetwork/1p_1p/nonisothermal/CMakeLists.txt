# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_md_boundary_ff1pni_pnm1pni
               SOURCES main.cc
               LABELS multidomain multidomain_boundary freeflowpnm 1p navierstokes porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=0
               COMMAND ${DUMUX_RUNTEST}
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               CMD_ARGS  --script fuzzy
               --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1pni_pnm1pni_ff-reference.vtu
                       ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1pni_pnm1pni_ff-00010.vtu
                       ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_ff1pni_pnm1pni_pnm-reference.vtp
                       ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1pni_pnm1pni_pnm-00010.vtp
               --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_ff1pni_pnm1pni paramsni.input
                          -Vtk.OutputName test_md_boundary_ff1pni_pnm1pni")
