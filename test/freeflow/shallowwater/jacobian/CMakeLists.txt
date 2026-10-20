# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES params.input)

# the Jacobian assembled with the two-ring stencil equals a finite-difference Jacobian of the global residual,
# up to round-off without limiter and up to the effect of the kinks of the limiter with it
dumux_add_test(NAME test_shallowwater_jacobian_secondorder
               SOURCES main.cc
               LABELS shallowwater unit
               CMD_ARGS params.input
                        -ShallowWater.Limiter None
                        -Problem.JacobianTolerance 1e-8)

dumux_add_test(NAME test_shallowwater_jacobian_secondorder_limited
               TARGET test_shallowwater_jacobian_secondorder
               LABELS shallowwater unit
               COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_jacobian_secondorder
               CMD_ARGS params.input
                        -ShallowWater.Limiter VanAlbada
                        -Problem.JacobianTolerance 1e-5)
