# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# Copy mesh files to build directory.
#  - sphere_eighth.msh         : moderate ~854-node mesh, the benchmark default
#                                (good accuracy vs the analytical solution).
#  - sphere_eighth_coarse.msh  : ~223-node mesh used only for the fast CI smoke test.
dune_symlink_to_source_files(FILES
    sphere_eighth.msh
    sphere_eighth_coarse.msh
    params.input
    plot_cryer.py
    run_benchmark.py
    make_setup_svg.py
)

# Build the Cryer benchmark executable. The CI test is a fast smoke test: it runs a
# few Implicit-Euler steps on the coarse mesh (the moderate default mesh + DIRK3 is
# for the actual accuracy comparison, run via plot_cryer.py).
dumux_add_test(NAME test_cryer_poroelastic_large_def
    SOURCES main.cc
    COMPILE_DEFINITIONS HAVE_UMFPACK=1
    CMAKE_GUARD "( dune-alugrid_FOUND AND SuiteSparse_UMFPACK_FOUND )"
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_cryer_poroelastic_large_def
    CMD_ARGS -Problem.PressureLoad 100.0
             -Grid.File sphere_eighth_coarse.msh
             -TimeLoop.Scheme ImplicitEuler
             -TimeLoop.TEnd 5.0
             -TimeLoop.Dt 1.0
             -TimeLoop.MaxTimeStepSize 1.0
             -SpatialParams.PermeabilityModel Constant
    TIMEOUT 300
    LABELS poromechanics)
