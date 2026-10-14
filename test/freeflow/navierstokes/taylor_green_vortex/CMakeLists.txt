# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES benchmark.py visualize.py animate.py verify.py)

# the momentum discretizations (the mass balance is always discretized with the Box scheme);
# pq2hybridhoquad is hybrid PQ2 with higher-order quadrature rules for the source term on the
# sub-control volumes and for the mass fluxes instead of the midpoint rule (for comparison)
set(TAYLORGREEN_SCHEMES_2 pq1bubble pq1bubblehybrid pq2hybrid pq2hybridhoquad pq3hybrid)
set(TAYLORGREEN_TYPETAG_pq1bubble TaylorGreenTestMomentumPQ1Bubble)
set(TAYLORGREEN_TYPETAG_pq1bubblehybrid TaylorGreenTestMomentumPQ1BubbleHybrid)
set(TAYLORGREEN_TYPETAG_pq2hybrid TaylorGreenTestMomentumPQ2Hybrid)
set(TAYLORGREEN_TYPETAG_pq2hybridhoquad TaylorGreenTestMomentumPQ2HybridHOQuad)
set(TAYLORGREEN_TYPETAG_pq3hybrid TaylorGreenTestMomentumPQ3Hybrid)
set(TAYLORGREEN_MASS_pq1bubble TaylorGreenTestMassBox)
set(TAYLORGREEN_MASS_pq1bubblehybrid TaylorGreenTestMassBox)
set(TAYLORGREEN_MASS_pq2hybrid TaylorGreenTestMassBox)
set(TAYLORGREEN_MASS_pq2hybridhoquad TaylorGreenTestMassBoxHOQuad)
set(TAYLORGREEN_MASS_pq3hybrid TaylorGreenTestMassBox)

# non-hybrid PQ1Bubble (one bubble per element) with Box pressure is only stable on simplices
set(TAYLORGREEN_SIMPLEX_pq1bubble 1)
set(TAYLORGREEN_SIMPLEX_pq1bubblehybrid 0)
set(TAYLORGREEN_SIMPLEX_pq2hybrid 0)
set(TAYLORGREEN_SIMPLEX_pq2hybridhoquad 0)
set(TAYLORGREEN_SIMPLEX_pq3hybrid 0)
set(TAYLORGREEN_GUARD_0 "HAVE_UMFPACK")
set(TAYLORGREEN_GUARD_1 "( HAVE_UMFPACK AND dune-alugrid_FOUND )")

# small grids for the regression tests (on coarser grids, the Newton solver does not
# converge for the stationary problem with PQ1Bubble)
set(TAYLORGREEN_CELLS_2 "16 16")
set(TAYLORGREEN_CELLS_2_pq1bubble "32 32")

foreach(dim 2)
  foreach(scheme ${TAYLORGREEN_SCHEMES_${dim}})
    set(name test_ff_navierstokes_taylorgreen_${dim}d_${scheme})
    set(simplex ${TAYLORGREEN_SIMPLEX_${scheme}})
    set(cells ${TAYLORGREEN_CELLS_${dim}})
    if(DEFINED TAYLORGREEN_CELLS_${dim}_${scheme})
      set(cells ${TAYLORGREEN_CELLS_${dim}_${scheme}})
    endif()

    # stationary problem (sustained by a source term)
    dumux_add_test(NAME ${name}
                   SOURCES main.cc
                   LABELS freeflow navierstokes taylorgreen
                   CMAKE_GUARD "${TAYLORGREEN_GUARD_${simplex}}"
                   TIMEOUT 900
                   COMPILE_DEFINITIONS TAYLORGREEN_DIM=${dim}
                   COMPILE_DEFINITIONS TYPETAG_MOMENTUM=${TAYLORGREEN_TYPETAG_${scheme}}
                   COMPILE_DEFINITIONS TYPETAG_MASS=${TAYLORGREEN_MASS_${scheme}}
                   COMPILE_DEFINITIONS SIMPLEX_GRID=${simplex}
                   COMMAND ${DUMUX_RUNTEST}
                   CMD_ARGS --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/${name}-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/${name}-00001.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/${name} params_${dim}d.input
                            -Problem.Name ${name}
                            -Problem.IsStationary true
                            -Grid.Cells '${cells}'")

    # instationary problem (decaying vortex), four time steps
    dumux_add_test(NAME ${name}_instationary
                   TARGET ${name}
                   LABELS freeflow navierstokes taylorgreen
                   CMAKE_GUARD "${TAYLORGREEN_GUARD_${simplex}}"
                   COMMAND ${DUMUX_RUNTEST}
                   CMD_ARGS --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/${name}_instationary-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/${name}_instationary-00004.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/${name} params_${dim}d.input
                            -Problem.Name ${name}_instationary
                            -Problem.IsStationary false
                            -Grid.Cells '${cells}'
                            -TimeLoop.TEnd 0.2")
  endforeach()
endforeach()

# spatial convergence tests (stationary problem)
foreach(scheme ${TAYLORGREEN_SCHEMES_2})
  dumux_add_test(NAME test_ff_navierstokes_taylorgreen_2d_${scheme}_convergence
                 TARGET test_ff_navierstokes_taylorgreen_2d_${scheme}
                 LABELS freeflow navierstokes taylorgreen
                 CMAKE_GUARD "${TAYLORGREEN_GUARD_${TAYLORGREEN_SIMPLEX_${scheme}}}"
                 TIMEOUT 1200
                 COMMAND ./benchmark.py
                 CMD_ARGS --dim 2 --schemes ${scheme} --study spatial --test)
endforeach()

# higher-order time integration (the instationary problem uses the multi-stage time stepping,
# selected with TimeLoop.Scheme: ImplicitEuler, CrankNicolson, DIRK2, DIRK3)
dumux_add_test(NAME test_ff_navierstokes_taylorgreen_2d_pq2hybrid_dirk3
               TARGET test_ff_navierstokes_taylorgreen_2d_pq2hybrid
               LABELS freeflow navierstokes taylorgreen
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${DUMUX_RUNTEST}
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_taylorgreen_2d_pq2hybrid_dirk3-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_taylorgreen_2d_pq2hybrid_dirk3-00004.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_taylorgreen_2d_pq2hybrid params_2d.input
                        -Problem.Name test_ff_navierstokes_taylorgreen_2d_pq2hybrid_dirk3
                        -Problem.IsStationary false
                        -Grid.Cells '${TAYLORGREEN_CELLS_2}'
                        -TimeLoop.Scheme DIRK3
                        -TimeLoop.TEnd 0.2")
