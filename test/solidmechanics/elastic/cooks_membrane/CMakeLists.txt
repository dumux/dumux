# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input cooksmembrane.msh)

foreach(_scheme box pq1bubble pq2 pq1fe pq2fe)
    if(_scheme STREQUAL "box")
        set(_model BoxModel)
    elseif(_scheme STREQUAL "pq1bubble")
        set(_model PQ1BubbleModel)
    elseif(_scheme STREQUAL "pq2")
        set(_model PQ2HybridModel)
    elseif(_scheme STREQUAL "pq1fe")
        set(_model PQ1FEModel)
    else()
        set(_model PQ2FEModel)
    endif()
    dumux_add_test(NAME test_cooks_membrane_${_scheme}
                   LABELS solidmechanics elastic
                   SOURCES main.cc
                   COMPILE_DEFINITIONS DISCRETIZATION_MODEL=${_model}
                   CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK AND DUMUX_HAVE_GRIDFORMAT )"
                   COMMAND ${DUMUX_RUNTEST}
                   CMD_ARGS --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_cooks_membrane_${_scheme}-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_${_scheme}.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_cooks_membrane_${_scheme} params.input -Problem.Name test_cooks_membrane_${_scheme}")
endforeach()
