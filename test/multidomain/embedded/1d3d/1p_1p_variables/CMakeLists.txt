# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

set(PREFIX test_md_embedded_1d3d_1p1p_variables)
set(RUN_ARGS params.input -Vtk.EnableVtkOutput false -Tissue.Grid.Cells "10 10 10" -Vessel.Grid.Cells 10)

function(add_variables_executable NAME BULKTAG LOWDIMTAG MODE)
  dumux_add_test(NAME ${PREFIX}_${NAME}
                 LABELS multidomain multidomain_embedded 1p
                 SOURCES main.cc
                 COMPILE_DEFINITIONS BULKTYPETAG=${BULKTAG} LOWDIMTYPETAG=${LOWDIMTAG} COUPLINGMODE=Embedded1d3dCouplingMode::${MODE}
                 COMMAND ./${PREFIX}_${NAME}
                 CMAKE_GUARD dune-foamgrid_FOUND
                 CMD_ARGS ${RUN_ARGS} -Problem.Name ${PREFIX}_${NAME})
endfunction()

set(BULK_TAGS Box PQ1BubbleHybrid PQ2Hybrid PQ1FE PQ2FE)
set(BULK_NAMES box pq1bubblehybrid pq2hybrid pq1fe pq2fe)
list(LENGTH BULK_TAGS NUM_BULK)
math(EXPR LAST_BULK "${NUM_BULK} - 1")
foreach(I RANGE ${LAST_BULK})
  list(GET BULK_TAGS ${I} TAG)
  list(GET BULK_NAMES ${I} BULK)
  add_variables_executable(${BULK}_box_average Tissue${TAG} BloodFlowBox Average)
  add_variables_executable(${BULK}_box_kernel Tissue${TAG} BloodFlowBox Kernel)
endforeach()
add_variables_executable(pq2hybrid_pq2hybrid_average TissuePQ2Hybrid BloodFlowPQ2Hybrid Average)
add_variables_executable(pq2hybrid_pq2hybrid_kernel TissuePQ2Hybrid BloodFlowPQ2Hybrid Kernel)
add_variables_executable(box_box_surface TissueBox BloodFlowBox Surface)
add_variables_executable(box_box_line TissueBox BloodFlowBox Line)
add_variables_executable(pq2hybrid_box_surface TissuePQ2Hybrid BloodFlowBox Surface)
add_variables_executable(pq2hybrid_box_line TissuePQ2Hybrid BloodFlowBox Line)

set(ALL_VARIANTS
    box_box_average box_box_kernel pq1bubblehybrid_box_average pq1bubblehybrid_box_kernel
    pq2hybrid_box_average pq2hybrid_box_kernel pq1fe_box_average pq1fe_box_kernel
    pq2fe_box_average pq2fe_box_kernel pq2hybrid_pq2hybrid_average pq2hybrid_pq2hybrid_kernel
    box_box_surface box_box_line pq2hybrid_box_surface pq2hybrid_box_line)

foreach(VARIANT ${ALL_VARIANTS})
  dumux_add_test(NAME ${PREFIX}_${VARIANT}_jacobian
                 LABELS multidomain multidomain_embedded 1p
                 TARGET ${PREFIX}_${VARIANT}
                 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PREFIX}_${VARIANT}
                 CMAKE_GUARD dune-foamgrid_FOUND
                 CMD_ARGS params.input -Vtk.EnableVtkOutput false -Test.CheckJacobian true
                          -Tissue.Grid.Cells "6 6 6" -Vessel.Grid.Cells 6 -MixedDimension.NumCircleSegments 16
                          -MixedDimension.KernelIntegrationCRL 0.1 -Problem.Name ${PREFIX}_${VARIANT}_jacobian)
endforeach()

foreach(VARIANT pq2hybrid_box_average pq1fe_box_kernel pq2hybrid_pq2hybrid_average)
  dumux_add_test(NAME ${PREFIX}_${VARIANT}_multithreading
                 LABELS multidomain multidomain_embedded 1p
                 THREADS 4
                 TARGET ${PREFIX}_${VARIANT}
                 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PREFIX}_${VARIANT}
                 CMAKE_GUARD dune-foamgrid_FOUND
                 CMD_ARGS params.input -Vtk.EnableVtkOutput false -Test.CheckMultithreading true -Test.CheckColoring true
                          -Tissue.Grid.Cells "8 8 8" -Vessel.Grid.Cells 8 -MixedDimension.NumCircleSegments 16
                          -MixedDimension.KernelIntegrationCRL 0.1 -Problem.Name ${PREFIX}_${VARIANT}_multithreading)
endforeach()

foreach(VARIANT pq1bubblehybrid_box_average pq1bubblehybrid_box_kernel pq2hybrid_box_average pq2hybrid_box_kernel
                pq1fe_box_average pq1fe_box_kernel pq2fe_box_average pq2fe_box_kernel
                pq2hybrid_pq2hybrid_average pq2hybrid_pq2hybrid_kernel)
  string(REGEX MATCH "(average|kernel)$" METHOD ${VARIANT})
  if(METHOD STREQUAL "kernel")
    set(RATE_ARGS --levels "10 20" "20 40" --min-rate 1.75 --tolerance 0.25)
  elseif(VARIANT MATCHES "^pq2hybrid")
    set(RATE_ARGS --levels "4 8" "8 16" "16 32" --tolerance 0.05)
  else()
    set(RATE_ARGS --levels "4 8" "8 16" "16 32" --tolerance 0.6 --max-error-ratio 3)
  endif()
  dumux_add_test(NAME ${PREFIX}_${VARIANT}_convergence
                 LABELS multidomain multidomain_embedded 1p
                 TARGET ${PREFIX}_${VARIANT}
                 COMMAND python3
                 CMAKE_GUARD dune-foamgrid_FOUND
                 CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../1p_1p/compare_rates.py
                          --command-a "${CMAKE_CURRENT_BINARY_DIR}/${PREFIX}_${VARIANT} params.input -Problem.Name ${PREFIX}_${VARIANT}_convergence -SpatialParams.Radius 0.1 -Vtk.EnableVtkOutput false"
                          --command-b "${CMAKE_CURRENT_BINARY_DIR}/${PREFIX}_box_box_${METHOD} params.input -Problem.Name ${PREFIX}_box_box_${METHOD}_convergence -SpatialParams.Radius 0.1 -Vtk.EnableVtkOutput false"
                          --errors-a ${PREFIX}_${VARIANT}_convergence-errors.csv
                          --errors-b ${PREFIX}_box_box_${METHOD}_convergence-errors.csv
                          --level-args "-Vessel.Grid.Cells {0} -Tissue.Grid.Cells '{1} {1} {0}'"
                          ${RATE_ARGS})
  add_dependencies(${PREFIX}_${VARIANT} ${PREFIX}_box_box_${METHOD})
endforeach()

foreach(METHOD average line surface)
  dumux_add_test(NAME ${PREFIX}_box_box_${METHOD}_equals_volume_variables
                 LABELS multidomain multidomain_embedded 1p
                 TARGET ${PREFIX}_box_box_${METHOD}
                 COMMAND python3
                 CMAKE_GUARD dune-foamgrid_FOUND
                 CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../1p_1p/compare_errors.py
                          --command-a "${CMAKE_CURRENT_BINARY_DIR}/${PREFIX}_box_box_${METHOD} params.input -Problem.Name ${PREFIX}_box_box_${METHOD}_equals -Vtk.EnableVtkOutput false -Problem.DirichletOnly true -SpatialParams.PermeabilityVariation 0"
                          --command-b "${CMAKE_CURRENT_BINARY_DIR}/../1p_1p/test_md_embedded_1d3d_1p1p_boxbox_${METHOD} params.input -Problem.Name ${PREFIX}_boxbox_${METHOD}_classic -Vtk.EnableVtkOutput false -Problem.DirichletOnly true -SpatialParams.PermeabilityVariation 0 -Test.DirectSolver true"
                          --errors-a ${PREFIX}_box_box_${METHOD}_equals-errors.csv
                          --errors-b ${PREFIX}_boxbox_${METHOD}_classic-errors.csv
                          --columns "L2Rel(p_3d)" "L2Rel(p_1d)"
                          --tolerance 1e-8)
  add_dependencies(${PREFIX}_box_box_${METHOD} test_md_embedded_1d3d_1p1p_boxbox_${METHOD})
endforeach()

dune_symlink_to_source_files(FILES "params.input")
