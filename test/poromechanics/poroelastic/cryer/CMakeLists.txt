# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

dumux_add_test(NAME test_cryer_poroelastic_small_def
    SOURCES main.cc
    CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_cryer_poroelastic_small_def
    CMD_ARGS params.input
             -Grid.File ../../sphere_eighth_coarse.msh
             -TimeLoop.Scheme ImplicitEuler
             -TimeLoop.TEnd 5.0
             -TimeLoop.Dt 1.0
             -TimeLoop.MaxTimeStepSize 1.0
             -Problem.CenterPressureTolerance 0.07
    TIMEOUT 300
    LABELS poromechanics)

dumux_add_test(NAME test_cryer_poroelastic_small_def_compressible
    TARGET test_cryer_poroelastic_small_def
    CMAKE_GUARD "( dune-alugrid_FOUND AND HAVE_UMFPACK )"
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_cryer_poroelastic_small_def
    CMD_ARGS params.input
             -Grid.File ../../sphere_eighth_coarse.msh
             -TimeLoop.Scheme ImplicitEuler
             -TimeLoop.TEnd 180.0
             -TimeLoop.Dt 2.5
             -TimeLoop.DtGrowthFactor 1.3
             -TimeLoop.MaxTimeStepSize 25.0
             -SpatialParams.BiotCoefficient 0.8
             -SpatialParams.StorageCoefficient 5e-7
             -Problem.CenterPressureTolerance 0.07
    TIMEOUT 300
    LABELS poromechanics)
