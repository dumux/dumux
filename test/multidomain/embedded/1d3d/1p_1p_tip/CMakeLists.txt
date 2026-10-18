# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input" "run_benchmark.py")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tip_tpfatpfa_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./test_md_embedded_1d3d_1p1p_tip_tpfatpfa_average
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tip_tpfatpfa_average)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tip_tpfatpfa_surface
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./test_md_embedded_1d3d_1p1p_tip_tpfatpfa_surface
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tip_tpfatpfa_surface)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tip_tpfatpfa_kernel
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./test_md_embedded_1d3d_1p1p_tip_tpfatpfa_kernel
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tip_tpfatpfa_kernel)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tip_boxtpfa_kernel
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./test_md_embedded_1d3d_1p1p_tip_boxtpfa_kernel
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tip_boxtpfa_kernel)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tip_convergence
              LABELS multidomain multidomain_embedded 1p
              TIMEOUT 1800
              TARGET test_md_embedded_1d3d_1p1p_tip_tpfatpfa_kernel
              COMMAND ./run_benchmark.py
              CMAKE_GUARD dune-foamgrid_FOUND)

set_tests_properties(test_md_embedded_1d3d_1p1p_tip_convergence PROPERTIES DEPENDS "test_md_embedded_1d3d_1p1p_tip_tpfatpfa_average;test_md_embedded_1d3d_1p1p_tip_tpfatpfa_surface;test_md_embedded_1d3d_1p1p_tip_boxtpfa_kernel")
