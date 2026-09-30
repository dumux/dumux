# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES benchmark.py)

# the momentum discretizations (the mass balance is always discretized with the Box scheme)
set(TAYLORGREEN_SCHEMES pq1bubble pq1bubblehybrid pq2hybrid)
set(TAYLORGREEN_TYPETAG_pq1bubble TaylorGreenTestMomentumPQ1Bubble)
set(TAYLORGREEN_TYPETAG_pq1bubblehybrid TaylorGreenTestMomentumPQ1BubbleHybrid)
set(TAYLORGREEN_TYPETAG_pq2hybrid TaylorGreenTestMomentumPQ2Hybrid)

# non-hybrid PQ1Bubble (one bubble per element) with Box pressure is only stable on simplices
set(TAYLORGREEN_SIMPLEX_pq1bubble 1)
set(TAYLORGREEN_SIMPLEX_pq1bubblehybrid 0)
set(TAYLORGREEN_SIMPLEX_pq2hybrid 0)
set(TAYLORGREEN_GUARD_0 "HAVE_UMFPACK")
set(TAYLORGREEN_GUARD_1 "( HAVE_UMFPACK AND dune-alugrid_FOUND )")

# small grids for the regression tests
set(TAYLORGREEN_CELLS_2 "16 16")

foreach(dim 2)
  foreach(scheme ${TAYLORGREEN_SCHEMES})
    set(name test_ff_navierstokes_taylorgreen_${dim}d_${scheme})
    set(simplex ${TAYLORGREEN_SIMPLEX_${scheme}})

    # stationary problem (sustained by a source term)
    dumux_add_test(NAME ${name}
                   SOURCES main.cc
                   LABELS freeflow navierstokes taylorgreen
                   CMAKE_GUARD "${TAYLORGREEN_GUARD_${simplex}}"
                   COMPILE_DEFINITIONS TAYLORGREEN_DIM=${dim}
                   COMPILE_DEFINITIONS TYPETAG_MOMENTUM=${TAYLORGREEN_TYPETAG_${scheme}}
                   COMPILE_DEFINITIONS SIMPLEX_GRID=${simplex}
                   COMMAND ${DUMUX_RUNTEST}
                   CMD_ARGS --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/${name}-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/${name}-00001.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/${name} params_${dim}d.input
                            -Problem.Name ${name}
                            -Problem.IsStationary true
                            -Grid.Cells '${TAYLORGREEN_CELLS_${dim}}'")

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
                            -Grid.Cells '${TAYLORGREEN_CELLS_${dim}}'
                            -TimeLoop.TEnd 0.2")
  endforeach()
endforeach()

# spatial convergence tests (stationary problem)
foreach(scheme ${TAYLORGREEN_SCHEMES})
  dumux_add_test(NAME test_ff_navierstokes_taylorgreen_2d_${scheme}_convergence
                 TARGET test_ff_navierstokes_taylorgreen_2d_${scheme}
                 LABELS freeflow navierstokes taylorgreen
                 CMAKE_GUARD "${TAYLORGREEN_GUARD_${TAYLORGREEN_SIMPLEX_${scheme}}}"
                 TIMEOUT 1200
                 COMMAND ./benchmark.py
                 CMD_ARGS --dim 2 --schemes ${scheme} --study spatial --test)
endforeach()

# temporal convergence test (implicit Euler)
dumux_add_test(NAME test_ff_navierstokes_taylorgreen_2d_pq2hybrid_temporal_convergence
               TARGET test_ff_navierstokes_taylorgreen_2d_pq2hybrid
               LABELS freeflow navierstokes taylorgreen
               CMAKE_GUARD HAVE_UMFPACK
               TIMEOUT 1200
               COMMAND ./benchmark.py
               CMD_ARGS --dim 2 --schemes pq2hybrid --study temporal --test)
