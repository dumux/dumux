# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input params_sweep.input)

dumux_add_test(NAME test_numericdifferentiation_complex
               SOURCES test_numericdifferentiation_complex.cc
               LABELS unit experimental)

dumux_add_test(NAME test_complex_helmholtz_box_bicgstab
               SOURCES main_manufactured.cc
               LABELS unit experimental
               COMPILE_DEFINITIONS TYPETAG=ComplexHelmholtzBox LINEARSOLVER=ILUBiCGSTABIstlSolver
               CMD_ARGS params.input)

dumux_add_test(NAME test_complex_helmholtz_box_umfpack
               SOURCES main_manufactured.cc
               LABELS unit experimental
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=ComplexHelmholtzBox LINEARSOLVER=UMFPackIstlSolver
               CMD_ARGS params.input)

# compiled with assertions enabled to cover the epsilon checks
dumux_add_test(NAME test_complex_helmholtz_box_umfpack_debug
               SOURCES main_manufactured.cc
               LABELS unit experimental
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=ComplexHelmholtzBox LINEARSOLVER=UMFPackIstlSolver
               COMPILE_FLAGS -UNDEBUG
               CMD_ARGS params.input)

dumux_add_test(NAME test_complex_helmholtz_tpfa_umfpack
               SOURCES main_manufactured.cc
               LABELS unit experimental
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=ComplexHelmholtzTpfa LINEARSOLVER=UMFPackIstlSolver
               CMD_ARGS params.input)

dumux_add_test(NAME test_complex_helmholtz_sweep
               SOURCES main_sweep.cc
               LABELS unit experimental
               CMD_ARGS params_sweep.input -Problem.MaxMode 3 -Grid.Cells "40 40")

add_subdirectory(facet)
