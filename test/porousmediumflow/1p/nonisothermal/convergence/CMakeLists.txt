# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py)

# the schemes with variables per local dof, each for conduction only (zero pressure gradient) and with convection;
# the hybrid PQ2 scheme converges with order 2 as for single-phase flow, and with convection the enthalpy
# upwinded between the two control volumes of a face limits all schemes to order 2
foreach(_scheme box pq1bubble pq2 pq3)
    if(_scheme STREQUAL "box")
        set(_typeTag OnePNIConvergenceBox)
        set(_grid "Dune::YaspGrid<2>")
        set(_rate 2)
        set(_guard HAVE_UMFPACK)
    elseif(_scheme STREQUAL "pq1bubble")
        set(_typeTag OnePNIConvergencePQ1Bubble)
        set(_grid "Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>")
        set(_rate 2)
        set(_guard "( HAVE_UMFPACK AND dune-alugrid_FOUND )")
    elseif(_scheme STREQUAL "pq2")
        set(_typeTag OnePNIConvergencePQ2)
        set(_grid "Dune::YaspGrid<2>")
        set(_rate 2)
        set(_guard HAVE_UMFPACK)
    else()
        set(_typeTag OnePNIConvergencePQ3)
        set(_grid "Dune::YaspGrid<2>")
        set(_rate 4)
        set(_guard HAVE_UMFPACK)
    endif()

    dumux_add_test(NAME test_1pni_convergence_${_scheme}
                   SOURCES main.cc
                   COMPILE_DEFINITIONS TYPETAG=${_typeTag} GRIDTYPE=${_grid}
                   CMAKE_GUARD ${_guard}
                   COMPILE_ONLY)

    dumux_add_test(NAME test_1pni_convergence_conduction_${_scheme}
                   TARGET test_1pni_convergence_${_scheme}
                   CMAKE_GUARD ${_guard}
                   LABELS porousmediumflow 1p 1pni
                   COMMAND ./convergencetest.py
                   CMD_ARGS test_1pni_convergence_${_scheme} test_1pni_convergence_conduction_${_scheme}
                            --expected-rate ${_rate}
                            params.input -Problem.PressureGradient "0.0 0.0")

    dumux_add_test(NAME test_1pni_convergence_convection_${_scheme}
                   TARGET test_1pni_convergence_${_scheme}
                   CMAKE_GUARD ${_guard}
                   LABELS porousmediumflow 1p 1pni
                   COMMAND ./convergencetest.py
                   CMD_ARGS test_1pni_convergence_${_scheme} test_1pni_convergence_convection_${_scheme}
                            --expected-rate 2
                            params.input -Problem.PressureGradient "-1.0 -0.5")
endforeach()
