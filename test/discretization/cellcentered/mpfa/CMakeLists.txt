# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_mpfafvgeometry
              SOURCES test_mpfafvgeometry.cc
              LABELS unit discretization)

dumux_add_test(NAME test_mpfafvgeometry_caching
              SOURCES test_mpfafvgeometry.cc
              LABELS unit discretization
              COMPILE_DEFINITIONS ENABLECACHING=1)
