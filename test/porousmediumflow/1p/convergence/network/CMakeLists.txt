# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES "convergencetest.py" "junction.dgf")

# with vertex control volumes on the median dual, the hybrid schemes of even order
# converge one order below the optimal rate in L2
dumux_add_test(NAME test_1p_convergence_line_pq2
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ2
               COMPILE_DEFINITIONS GRIDTYPE=Dune::OneDGrid
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_line_pq2 line_pq2
                        --l2-rate 2 --h1-rate 2 params_line.input)

dumux_add_test(NAME test_1p_convergence_line_pq3
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ3
               COMPILE_DEFINITIONS GRIDTYPE=Dune::OneDGrid
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_line_pq3 line_pq3
                        --l2-rate 4 --h1-rate 3 params_line.input)

dumux_add_test(NAME test_1p_convergence_junction_pq2
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ2
               COMPILE_DEFINITIONS GRIDTYPE=Dune::FoamGrid<1,3>
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_junction_pq2 junction_pq2
                        --l2-rate 2 --h1-rate 2 params_junction.input)

dumux_add_test(NAME test_1p_convergence_junction_pq3
               SOURCES main.cc
               CMAKE_GUARD "( HAVE_UMFPACK AND dune-foamgrid_FOUND )"
               COMPILE_DEFINITIONS TYPETAG=OnePNetworkConvergencePQ3
               COMPILE_DEFINITIONS GRIDTYPE=Dune::FoamGrid<1,3>
               LABELS porousmediumflow 1p
               COMMAND ./convergencetest.py
               CMD_ARGS test_1p_convergence_junction_pq3 junction_pq3
                        --l2-rate 4 --h1-rate 3 params_junction.input)
