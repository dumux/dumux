# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_boxfvgeometry
              SOURCES test_boxfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false
              LABELS unit discretization)

dumux_add_test(NAME test_boxfvgeometry_caching
              SOURCES test_boxfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true
              LABELS unit discretization)
