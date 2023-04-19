# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_diamondgridgeometry
               SOURCES test_diamondgridgeometry.cc
               COMPILE_DEFINITIONS ENABLE_CACHING=true
               LABELS unit discretization)
