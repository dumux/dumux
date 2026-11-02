# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# The momentum balance with a prescribed density that is transported by the same stencil as the
# momentum is solved up to the discretization error of the hybrid dofs, for the form of the advective flux
# that matches the stencil.
dumux_add_test(NAME test_ff_navierstokes_density_advection
              SOURCES main.cc
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_density_advection
              CMD_ARGS params.input
                       -Study.DensityField upwind
                       -Study.DensityRatios "10 1000"
                       -Study.MaxError "1e-3 1e-2"
                       -Study.MinOrder "2.0 2.0")

dumux_add_test(NAME test_ff_navierstokes_density_advection_interpolated_density
              TARGET test_ff_navierstokes_density_advection
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_density_advection
              CMD_ARGS params.input
                       -Flux.UpwindMomentumDensity false
                       -Study.DensityField central
                       -Study.DensityRatios "10"
                       -Study.MaxError "2e-4"
                       -Study.MinOrder "2.5")

dumux_add_test(NAME test_ff_navierstokes_density_advection_constant_density
              TARGET test_ff_navierstokes_density_advection
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_density_advection
              CMD_ARGS params.input
                       -Study.BaseDensity 1000
                       -Study.NumLevels 2
                       -Study.DensityRatios "1"
                       -Study.MaxError "1e-12")

dumux_add_test(NAME test_ff_navierstokes_density_advection_constant_density_interpolated
              TARGET test_ff_navierstokes_density_advection
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_density_advection
              CMD_ARGS params.input
                       -Flux.UpwindMomentumDensity false
                       -Study.BaseDensity 1000
                       -Study.NumLevels 2
                       -Study.DensityRatios "1"
                       -Study.MaxError "1e-12")
