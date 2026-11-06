# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input grids)

# meshes of the convergence tests with the given number of elements per side
if(Gmsh_FOUND)
  foreach(cells 10 20 30 40 50 60 70 80)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hybridgrid_${cells}.msh
                       COMMAND ${GMSH_EXECUTABLE} -2 -format msh2 -setnumber numElemsPerSide ${cells} ${CMAKE_CURRENT_SOURCE_DIR}/grids/hybridgrid.geo -o ${CMAKE_CURRENT_BINARY_DIR}/hybridgrid_${cells}.msh
                       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/grids/hybridgrid.geo)
    list(APPEND COMPLEX_FACET_MESHES ${CMAKE_CURRENT_BINARY_DIR}/hybridgrid_${cells}.msh)
  endforeach()
  add_custom_target(test_complex_facet_helmholtz_meshes DEPENDS ${COMPLEX_FACET_MESHES})
endif()

# complex amplitude and complex squared wave number
dumux_add_test(NAME test_complex_facet_helmholtz
               SOURCES main.cc
               LABELS multidomain multidomain_facet experimental
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND Gmsh_FOUND AND HAVE_UMFPACK )"
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_complex_facet_helmholtz params.input -Problem.Name test_complex_facet_helmholtz -Problem.Amplitude '0.6 -0.8' -Problem.WaveNumberSquared '20 -5'"
                        --levels 10 20 30 40 50 60 70 80
                        --level-args "-Grid.File hybridgrid_{0}.msh -Grid.NumElemsPerSide {0}"
                        --errors test_complex_facet_helmholtz-errors.csv
                        --check "L2Rel(p_m) finest 1.9"
                        --check "L2Rel(p_f) finest 1.9")

# real data: the discrete problem of the analytical single-phase facet test
dumux_add_test(NAME test_complex_facet_helmholtz_real
               TARGET test_complex_facet_helmholtz
               LABELS multidomain multidomain_facet experimental
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND Gmsh_FOUND AND HAVE_UMFPACK )"
               COMMAND ${DUMUX_CONVERGENCETEST}
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_complex_facet_helmholtz params.input -Problem.Name test_complex_facet_helmholtz_real"
                        --levels 10 20 30 40 50 60 70 80
                        --level-args "-Grid.File hybridgrid_{0}.msh -Grid.NumElemsPerSide {0}"
                        --errors test_complex_facet_helmholtz_real-errors.csv
                        --check "L2Rel(p_m) finest 1.9"
                        --check "L2Rel(p_f) finest 1.9")

if(Gmsh_FOUND)
  add_dependencies(test_complex_facet_helmholtz test_complex_facet_helmholtz_meshes)
endif()
