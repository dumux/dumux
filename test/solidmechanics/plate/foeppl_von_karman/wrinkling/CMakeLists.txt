# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input wrinkling_benchmark.py plot_wrinkles.py)

dumux_add_test(NAME test_fvk_wrinkling
              LABELS solidmechanics plate
              SOURCES main.cc
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_fvk_wrinkling
              CMD_ARGS params.input
                       -Grid.Cells "48 16"
                       -Problem.LoadSteps 4096
                       -Relaxation.MaxSteps 150
                       -Problem.ShearResultantTolerance 0.05
                       -Problem.NormalAngleTolerance 10
                       -Problem.WavelengthTolerance 0.3
                       -Problem.MinAmplitudeOverThickness 2
                       -Problem.Name test_fvk_wrinkling)

# the sheet tensioned across its width and not sheared, against the uniform strip solution
dumux_add_test(NAME test_fvk_wrinkling_prestrain_flat
              TARGET test_fvk_wrinkling
              LABELS solidmechanics plate
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_fvk_wrinkling
              CMD_ARGS params.input
                       -Grid.Cells "96 32"
                       -Problem.LoadSteps 0
                       -Problem.PreStrain 1.0486e4
                       -Problem.SeedTensionField false
                       -Relaxation.MaxSteps 0
                       -Relaxation.FinalStaticSolve true
                       -Problem.PreStrainResultantTolerance 0.03
                       -Problem.Name test_fvk_wrinkling_prestrain_flat)

# sheared with a tenth of the shear strain locked in across the sheet, which turns the
# crests towards y by 2.9 degrees and raises the tension-field shear resultant by 10 %
dumux_add_test(NAME test_fvk_wrinkling_prestrain
              TARGET test_fvk_wrinkling
              LABELS solidmechanics plate
              CMAKE_GUARD dune-uggrid_FOUND
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_fvk_wrinkling
              CMD_ARGS params.input
                       -Grid.Cells "48 16"
                       -Problem.LoadSteps 4096
                       -Problem.PreStrain 409.6
                       -Relaxation.MaxSteps 150
                       -Problem.ShearResultantTolerance 0.05
                       -Problem.NormalAngleTolerance 2
                       -Problem.MinAmplitudeOverThickness 2
                       -Problem.Name test_fvk_wrinkling_prestrain)
