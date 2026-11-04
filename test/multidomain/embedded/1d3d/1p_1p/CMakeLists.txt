# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_average
              THREADS 4
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${DUMUX_RUNTEST}
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_average")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_boxtpfa_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${DUMUX_RUNTEST}
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxtpfa_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxtpfa_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxtpfa_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_boxtpfa_average -Tissue.Grid.Cells \"19 19 19\"")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfabox_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${DUMUX_RUNTEST}
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfabox_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_tpfatpfa_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfabox_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfabox_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfabox_average")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_boxbox_average
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=BloodFlowBox COUPLINGMODE=Embedded1d3dCouplingMode::Average
              COMMAND ${DUMUX_RUNTEST}
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_1d-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxbox_average_1d-00001.vtp
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_embedded_1d3d_1p1p_boxbox_average_3d-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxbox_average_3d-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxbox_average params.input \
                                   -Vtk.OutputName test_md_embedded_1d3d_1p1p_boxbox_average -Tissue.Grid.Cells \"19 19 19\"")

# make sure these configurations compile and run too
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_surface
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Surface
              COMMAND ./test_md_embedded_1d3d_1p1p_tpfatpfa_surface
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_surface)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_line
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Line
              COMMAND ./test_md_embedded_1d3d_1p1p_tpfatpfa_line
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_line)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_kernel
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              COMMAND ./test_md_embedded_1d3d_1p1p_tpfatpfa_kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfatpfa_kernel)

# The kernel (distributed source) method with vertex-centred dofs. The bulk source is spread
# over the sub-control volumes of the element containing each sample point, so the exchange
# term is a load vector rather than a single cell value.
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_boxtpfa_kernel
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              COMMAND ./test_md_embedded_1d3d_1p1p_boxtpfa_kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_boxtpfa_kernel)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_boxbox_kernel
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueBox LOWDIMTYPETAG=BloodFlowBox COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              COMMAND ./test_md_embedded_1d3d_1p1p_boxbox_kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_boxbox_kernel)

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfabox_kernel
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueCC LOWDIMTYPETAG=BloodFlowBox COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              COMMAND ./test_md_embedded_1d3d_1p1p_tpfabox_kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              CMD_ARGS params.input -Vtk.OutputName test_md_embedded_1d3d_1p1p_tpfabox_kernel)

# Unstructured bulk grid, which exercises the bounding-box-tree point location together with
# the sub-control-volume distribution. Box only: tpfa is not consistent on simplices, so a
# tpfa variant here would measure the scheme's own inconsistency rather than the coupling.
# The meshes are generated from box.geo at build time, so gmsh is required as well.
if(Gmsh_FOUND)
  foreach(clmax 0.30 0.20 0.15)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/box_${clmax}.msh
                       COMMAND ${GMSH_EXECUTABLE} -3 -format msh2 -clmax ${clmax} ${CMAKE_CURRENT_SOURCE_DIR}/box.geo -o ${CMAKE_CURRENT_BINARY_DIR}/box_${clmax}.msh
                       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/box.geo)
    list(APPEND EMBEDDED_BOX_MESHES ${CMAKE_CURRENT_BINARY_DIR}/box_${clmax}.msh)
  endforeach()
  add_custom_target(test_md_embedded_1d3d_1p1p_box_meshes DEPENDS ${EMBEDDED_BOX_MESHES})
endif()

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_boxtpfa_kernel_unstructured
              LABELS multidomain multidomain_embedded 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=TissueALUBox LOWDIMTYPETAG=BloodFlowCC COUPLINGMODE=Embedded1d3dCouplingMode::Kernel
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND Gmsh_FOUND )"
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_boxtpfa_kernel_unstructured params_unstructured.input -Problem.Name test_md_embedded_1d3d_1p1p_boxtpfa_kernel_unstructured"
                       --levels 0.30 0.20 0.15
                       --level-args "-Tissue.Grid.File box_{}.msh"
                       --errors test_md_embedded_1d3d_1p1p_boxtpfa_kernel_unstructured-errors.csv
                       --check "L2Rel(p_3d) finest 1.5")

if(Gmsh_FOUND)
  add_dependencies(test_md_embedded_1d3d_1p1p_boxtpfa_kernel_unstructured test_md_embedded_1d3d_1p1p_box_meshes)
endif()

# mean rates on 20, 40 and 80 tissue cells per direction within 0.05 of the rates
# of the convergence study reproduced by plot_convergence.py
dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_average_convergence
              LABELS multidomain multidomain_embedded 1p
              TIMEOUT 3600
              TARGET test_md_embedded_1d3d_1p1p_tpfatpfa_average
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_average params.input -Problem.Name test_md_embedded_1d3d_1p1p_tpfatpfa_average_convergence -SpatialParams.Radius 0.1 -Vtk.EnableVtkOutput false"
                       --levels "10 20" "20 40" "40 80"
                       --level-args "-Vessel.Grid.Cells {0} -Tissue.Grid.Cells '{1} {1} {0}'"
                       --errors test_md_embedded_1d3d_1p1p_tpfatpfa_average_convergence-errors.csv
                       --check "L2Rel(p_3d) mean 1.2983 1.3983"
                       --check "L2Rel(p_1d) mean 1.3045 1.4045"
                       --check "L2Rel(q) mean 1.7949 1.8949")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_surface_convergence
              LABELS multidomain multidomain_embedded 1p
              TIMEOUT 3600
              TARGET test_md_embedded_1d3d_1p1p_tpfatpfa_surface
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_surface params.input -Problem.Name test_md_embedded_1d3d_1p1p_tpfatpfa_surface_convergence -SpatialParams.Radius 0.1 -Vtk.EnableVtkOutput false"
                       --levels "10 20" "20 40" "40 80"
                       --level-args "-Vessel.Grid.Cells {0} -Tissue.Grid.Cells '{1} {1} {0}'"
                       --errors test_md_embedded_1d3d_1p1p_tpfatpfa_surface_convergence-errors.csv
                       --check "L2Rel(p_3d) mean 1.4193 1.5193"
                       --check "L2Rel(p_1d) mean 1.7300 1.8300"
                       --check "L2Rel(q) mean 2.2215 2.3215")

dumux_add_test(NAME test_md_embedded_1d3d_1p1p_tpfatpfa_kernel_convergence
              LABELS multidomain multidomain_embedded 1p
              TIMEOUT 3600
              TARGET test_md_embedded_1d3d_1p1p_tpfatpfa_kernel
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_embedded_1d3d_1p1p_tpfatpfa_kernel params.input -Problem.Name test_md_embedded_1d3d_1p1p_tpfatpfa_kernel_convergence -SpatialParams.Radius 0.1 -Vtk.EnableVtkOutput false"
                       --levels "10 20" "20 40" "40 80"
                       --level-args "-Vessel.Grid.Cells {0} -Tissue.Grid.Cells '{1} {1} {0}'"
                       --errors test_md_embedded_1d3d_1p1p_tpfatpfa_kernel_convergence-errors.csv
                       --check "L2Rel(p_3d) mean 2.0083 2.1083"
                       --check "L2Rel(p_1d) mean 2.0427 2.1427"
                       --check "L2Rel(q) mean 2.5322 2.6322")

dune_symlink_to_source_files(FILES "params.input" "params_unstructured.input"
                                   "plot_convergence.py" "box.geo")
