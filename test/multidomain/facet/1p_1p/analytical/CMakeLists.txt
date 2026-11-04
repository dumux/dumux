# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "grids" "params.input")

# meshes of the convergence tests with the given number of elements per side
if(Gmsh_FOUND)
  foreach(cells 10 20 30 40 50 60 70 80 90 100)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hybridgrid_${cells}.msh
                       COMMAND ${GMSH_EXECUTABLE} -2 -format msh2 -setnumber numElemsPerSide ${cells} ${CMAKE_CURRENT_SOURCE_DIR}/grids/hybridgrid.geo -o ${CMAKE_CURRENT_BINARY_DIR}/hybridgrid_${cells}.msh
                       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/grids/hybridgrid.geo)
    list(APPEND FACET_ANALYTICAL_MESHES ${CMAKE_CURRENT_BINARY_DIR}/hybridgrid_${cells}.msh)
  endforeach()
  add_custom_target(test_md_facet_1p1p_analytical_meshes DEPENDS ${FACET_ANALYTICAL_MESHES})
endif()

dumux_add_test(NAME test_md_facet_1p1p_tpfa_convergence
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND Gmsh_FOUND )"
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkTpfa LOWDIMTYPETAG=OnePLowDimTpfa
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_tpfa_convergence params.input -Problem.Name test_md_facet_1p1p_tpfa_convergence -Vtk.OutputName test_md_facet_1p1p_tpfa_convergence -LowDim.SpatialParams.Permeability 1e-4"
                       --levels 10 20 30 40 50 60 70 80 90 100
                       --level-args "-Grid.File hybridgrid_{0}.msh -Grid.NumElemsPerSide {0}"
                       --errors test_md_facet_1p1p_tpfa_convergence-errors.csv
                       --check "L2Rel(p_m) finest 1.9")

if(Gmsh_FOUND)
  add_dependencies(test_md_facet_1p1p_tpfa_convergence test_md_facet_1p1p_analytical_meshes)
endif()

dumux_add_test(NAME test_md_facet_1p1p_box_convergence
              LABELS multidomain multidomain_facet 1p
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND Gmsh_FOUND )"
              SOURCES main.cc
              COMPILE_DEFINITIONS BULKTYPETAG=OnePBulkBox LOWDIMTYPETAG=OnePLowDimBox
              COMMAND ${DUMUX_CONVERGENCETEST}
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_box_convergence params.input -Problem.Name test_md_facet_1p1p_box_convergence -Vtk.OutputName test_md_facet_1p1p_box_convergence -LowDim.SpatialParams.Permeability 1e-4"
                       --levels 10 20 30 40 50 60 70 80 90 100
                       --level-args "-Grid.File hybridgrid_{0}.msh -Grid.NumElemsPerSide {0}"
                       --errors test_md_facet_1p1p_box_convergence-errors.csv
                       --check "L2Rel(p_m) finest 1.9")

if(Gmsh_FOUND)
  add_dependencies(test_md_facet_1p1p_box_convergence test_md_facet_1p1p_analytical_meshes)
endif()
