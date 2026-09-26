# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input grids convergencetest.py)

# complex amplitude and complex squared wave number
dumux_add_test(NAME test_complex_facet_helmholtz
               SOURCES main.cc
               LABELS multidomain multidomain_facet experimental
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND Gmsh_FOUND AND HAVE_UMFPACK )"
               COMMAND ./convergencetest.py
               CMD_ARGS test_complex_facet_helmholtz test_complex_facet_helmholtz 1.9 -Problem.Amplitude "0.6 -0.8" -Problem.WaveNumberSquared "20 -5")

# real data: the discrete problem of the analytical single-phase facet test
dumux_add_test(NAME test_complex_facet_helmholtz_real
               TARGET test_complex_facet_helmholtz
               LABELS multidomain multidomain_facet experimental
               CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND AND Gmsh_FOUND AND HAVE_UMFPACK )"
               COMMAND ./convergencetest.py
               CMD_ARGS test_complex_facet_helmholtz test_complex_facet_helmholtz_real 1.9)
