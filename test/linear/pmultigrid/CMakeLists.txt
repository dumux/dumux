# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params_2d.input params_3d.input)

# add_pmultigrid_laplace_test(<name> <model> <grid> <solver> <params> <max iterations> <min rate> [MPI_RANKS ...] [GUARD <guard>])
function(add_pmultigrid_laplace_test name model grid solver params maxIt minRate)
  cmake_parse_arguments(ARG "" "GUARD" "MPI_RANKS" ${ARGN})
  set(guard TRUE)
  if(ARG_GUARD)
    set(guard ${ARG_GUARD})
  endif()
  set(ranks 1)
  if(ARG_MPI_RANKS)
    set(ranks ${ARG_MPI_RANKS})
  endif()
  dumux_add_test(NAME test_pmultigrid_laplace_${name}
                 SOURCES main.cc
                 COMPILE_DEFINITIONS SCHEME_MODEL=${model} GRIDTYPE=${grid} SOLVER=${solver}
                 CMAKE_GUARD ${guard}
                 LABELS linear unit
                 MPI_RANKS ${ranks}
                 TIMEOUT 600
                 CMD_ARGS ${params} -Test.MaxIterations ${maxIt} -Test.MinRate ${minRate})
endfunction()

set(simplex2d "Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>")
set(simplex3d "Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>")

add_pmultigrid_laplace_test(pq2fe_2d_simplex PQ2FEModel ${simplex2d} PMultigridCGIstlSolver params_2d.input 30 2.8
                            MPI_RANKS 1 2 4 GUARD dune-alugrid_FOUND)
add_pmultigrid_laplace_test(pq2fe_2d_simplex_gpcg PQ2FEModel ${simplex2d} PMultigridGeneralizedPCGIstlSolver params_2d.input 30 2.8
                            MPI_RANKS 1 2 4 GUARD dune-alugrid_FOUND)
add_pmultigrid_laplace_test(pq3fe_2d_simplex PQ3FEModel ${simplex2d} PMultigridCGIstlSolver params_2d.input 30 3.8
                            GUARD dune-alugrid_FOUND)
add_pmultigrid_laplace_test(pq2fe_3d_simplex PQ2FEModel ${simplex3d} PMultigridCGIstlSolver params_3d.input 30 2.8
                            MPI_RANKS 1 2 4 GUARD dune-alugrid_FOUND)
add_pmultigrid_laplace_test(pq2fe_2d_cube PQ2FEModel "Dune::YaspGrid<2>" PMultigridCGIstlSolver params_2d.input 30 2.8
                            MPI_RANKS 1 2 4)
add_pmultigrid_laplace_test(pq3fe_2d_cube PQ3FEModel "Dune::YaspGrid<2>" PMultigridCGIstlSolver params_2d.input 30 3.8)
add_pmultigrid_laplace_test(pq2fe_3d_cube PQ2FEModel "Dune::YaspGrid<3>" PMultigridCGIstlSolver params_3d.input 30 2.8
                            MPI_RANKS 1 4)
add_pmultigrid_laplace_test(pq2hybrid_2d_simplex PQ2HybridModel ${simplex2d} PMultigridBiCGSTABIstlSolver params_2d.input 30 1.8
                            MPI_RANKS 1 2 4 GUARD dune-alugrid_FOUND)


# overlapping decompositions (YaspGrid with overlap)
dumux_add_test(NAME test_pmultigrid_laplace_pq2fe_2d_cube_overlap
               TARGET test_pmultigrid_laplace_pq2fe_2d_cube
               LABELS linear unit
               MPI_RANKS 2 4
               TIMEOUT 600
               CMD_ARGS params_2d.input -Grid.Overlap 1 -Test.MaxIterations 30 -Test.MinRate 2.8)

dumux_add_test(NAME test_pmultigrid_laplace_pq2fe_3d_cube_overlap
               TARGET test_pmultigrid_laplace_pq2fe_3d_cube
               LABELS linear unit
               MPI_RANKS 4
               TIMEOUT 600
               CMD_ARGS params_3d.input -Grid.Overlap 1 -Test.MaxIterations 30 -Test.MinRate 2.8)
