# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input convergencetest.py)

# first-order flux: first order on the manufactured solution
dumux_add_test(NAME test_shallowwater_manufactured
               SOURCES main.cc
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=ManufacturedSolution LINEARSOLVER=UMFPackIstlSolver
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_shallowwater_manufactured 0.9)

# second-order flux: second order, with and without limiter, and a Newton step of order at least 1.8
# with the Jacobian of the two-ring stencil
dumux_add_test(NAME test_shallowwater_manufactured_secondorder
               SOURCES main.cc
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=ManufacturedSolutionSecondOrder LINEARSOLVER=UMFPackIstlSolver
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_shallowwater_manufactured_secondorder 1.9 --newton-order -ShallowWater.Limiter None)

# with the limiter, the steady state is reached by implicit Euler steps of growing size
dumux_add_test(NAME test_shallowwater_manufactured_secondorder_limited
               TARGET test_shallowwater_manufactured_secondorder
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               TIMEOUT 600
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_shallowwater_manufactured_secondorder 1.9 -ShallowWater.Limiter VanAlbada
                        -TimeLoop.TEnd 20 -TimeLoop.DtInitial 0.01 -TimeLoop.MaxTimeStepSize 2
                        -Newton.MaxRelativeShift 1e-10)

dumux_add_test(NAME test_shallowwater_manufactured_secondorder_simplex
               SOURCES main.cc
               LABELS shallowwater
               CMAKE_GUARD "( dune-uggrid_FOUND AND HAVE_UMFPACK )"
               COMPILE_DEFINITIONS GRIDTYPE=Dune::UGGrid<2> TYPETAG=ManufacturedSolutionSecondOrder LINEARSOLVER=UMFPackIstlSolver
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_shallowwater_manufactured_secondorder_simplex 1.9 -Grid.CellType Simplex -ShallowWater.Limiter None)

# time-dependent solution with BDF2: second order in space and time, with the time step size refined with the cell size
dumux_add_test(NAME test_shallowwater_manufactured_secondorder_bdf2
               SOURCES main_transient.cc
               LABELS shallowwater
               CMAKE_GUARD HAVE_UMFPACK
               COMPILE_DEFINITIONS TYPETAG=ManufacturedSolutionSecondOrder LINEARSOLVER=UMFPackIstlSolver
               TIMEOUT 1200
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/convergencetest.py
               CMD_ARGS test_shallowwater_manufactured_secondorder_bdf2 1.9 --dt-times-cells 1
                        -Problem.Frequency 2 -TimeLoop.TEnd 1 -TimeStepper.BDFOrder 2 -ShallowWater.Limiter None
                        -Newton.MaxRelativeShift 1e-10)
