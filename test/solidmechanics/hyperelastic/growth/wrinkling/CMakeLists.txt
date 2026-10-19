# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input wrinklingtest.py branchtest.py reference.py
                                   run_benchmark.py)

# critical growth stretch and number of half wavelengths from the linearized bilayer problem
set(_reference_args 1.1430689262 2)

dumux_add_test(NAME test_hyperelastic_growth_wrinkling_box
               LABELS solidmechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=FilmSubstrateBox
               CMAKE_GUARD "( DUMUX_HAVE_SPECTRA AND HAVE_UMFPACK )"
               COMMAND ./wrinklingtest.py
               CMD_ARGS test_hyperelastic_growth_wrinkling_box params.input ${_reference_args} 1.8)

dumux_add_test(NAME test_hyperelastic_growth_wrinkling_pq2hybrid
               LABELS solidmechanics
               SOURCES main.cc
               COMPILE_DEFINITIONS TYPETAG=FilmSubstratePQ2Hybrid
               CMAKE_GUARD "( DUMUX_HAVE_SPECTRA AND HAVE_UMFPACK )"
               COMMAND ./wrinklingtest.py
               CMD_ARGS test_hyperelastic_growth_wrinkling_pq2hybrid params.input ${_reference_args} 1.8)

dumux_add_test(NAME test_hyperelastic_growth_wrinkling_branch_box
               LABELS solidmechanics
               TARGET test_hyperelastic_growth_wrinkling_box
               CMAKE_GUARD "( DUMUX_HAVE_SPECTRA AND HAVE_UMFPACK )"
               COMMAND ./branchtest.py
               CMD_ARGS test_hyperelastic_growth_wrinkling_box params.input "0.01 0.02 0.04 0.08 0.16" 1e-6)
