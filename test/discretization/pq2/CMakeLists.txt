# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

dumux_add_test(NAME test_pq2fvgeometry
               SOURCES test_pq2fvgeometry.cc
               COMPILE_DEFINITIONS ENABLE_CACHING=true
               LABELS unit discretization)

set(CMAKE_BUILD_TYPE Debug)