# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# explicit tracer tests (mass fractions)
dumux_add_test(NAME test_tracer_explicit_tpfa
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa params.input -Problem.Name test_tracer_explicit_tpfa")

dumux_add_test(NAME test_tracer_explicit_box
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box params.input -Problem.Name test_tracer_explicit_box")

# explicit tracer tests (mole fractions, should yield same result)
dumux_add_test(NAME test_tracer_explicit_tpfa_mol
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa_mol-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_tpfa_mol params.input -Problem.Name test_tracer_explicit_tpfa_mol")

dumux_add_test(NAME test_tracer_explicit_box_mol
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_explicit_box-reference.vtu
                           ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box_mol-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_explicit_box_mol params.input -Problem.Name test_tracer_explicit_box_mol")

# implicit tracer tests
dumux_add_test(NAME test_tracer_implicit_tpfa
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=true USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_tpfa-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_tpfa params.input -Problem.Name test_tracer_implicit_tpfa")

dumux_add_test(NAME test_tracer_implicit_box
              LABELS porousmediumflow tracer
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=true USEMOLES=false
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_box params.input -Problem.Name test_tracer_implicit_box")

dumux_add_test(NAME test_tracer_implicit_dispersion_box
              LABELS porousmediumflow tracer dispersion
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=true USEMOLES=false ENABLEDISPERSION=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_dispersion_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_dispersion_box-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_dispersion_box params.input -Problem.Name test_tracer_implicit_dispersion_box
                       -Problem.AlphaL 0.02
                       -Problem.AlphaT 0.008")

dumux_add_test(NAME test_tracer_implicit_dispersion_tpfa
              LABELS porousmediumflow tracer dispersion
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=true USEMOLES=false ENABLEDISPERSION=true
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_dispersion_tpfa-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_dispersion_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_dispersion_tpfa params.input -Problem.Name test_tracer_implicit_dispersion_tpfa
                       -Problem.AlphaL 0.02
                       -Problem.AlphaT 0.008")

# implicit tracer tests with the direct solver MUMPS, which reuses the factorization of the constant
# matrix, sequentially and on two processes (box also on a non-overlapping decomposition); the
# velocity output is not computed in the parallel runs (not supported for cell-centered schemes
# on overlapping decompositions, averaged differently at the process border for box)
foreach(disc tpfa box)
    if(disc STREQUAL "tpfa")
        set(typetag TracerTestTpfa)
        set(overlaps 1)
    else()
        set(typetag TracerTestBox)
        set(overlaps 0 1)
    endif()
    dumux_add_test(NAME test_tracer_implicit_${disc}_mumps
                  LABELS porousmediumflow tracer
                  SOURCES main.cc
                  COMPILE_DEFINITIONS TYPETAG=${typetag} IMPLICIT=true USEMOLES=false USE_MUMPS_SOLVER=1
                  CMAKE_GUARD MUMPS_FOUND
                  COMMAND ${DUMUX_RUNTEST}
                  CMD_ARGS --script fuzzy
                           --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_${disc}-reference.vtu
                                   ${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_${disc}_mumps-00010.vtu
                           --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_${disc}_mumps params.input -Problem.Name test_tracer_implicit_${disc}_mumps")
    foreach(overlap ${overlaps})
        dumux_add_test(NAME test_tracer_implicit_${disc}_mumps_parallel_ovlp_${overlap}
                      TARGET test_tracer_implicit_${disc}_mumps
                      LABELS porousmediumflow tracer parallel
                      CMAKE_GUARD MUMPS_FOUND
                      MPI_RANKS 2
                      TIMEOUT 1200
                      COMMAND ${DUMUX_RUNTEST}
                      CMD_ARGS --script fuzzy --ignore "process rank" "velocity_Groundwater (m/s)"
                               --files ${CMAKE_SOURCE_DIR}/test/references/test_tracer_implicit_${disc}-reference.vtu
                                       ${CMAKE_CURRENT_BINARY_DIR}/s0002-test_tracer_implicit_${disc}_mumps_parallel_ovlp_${overlap}-00010.pvtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_implicit_${disc}_mumps params.input
                               -Problem.Name test_tracer_implicit_${disc}_mumps_parallel_ovlp_${overlap} -Grid.Overlap ${overlap}
                               -Vtk.AddVelocity false")
    endforeach()
endforeach()
