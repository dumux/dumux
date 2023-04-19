# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_pq1bubblefvgeometry
               SOURCES test_pq1bubblefvgeometry.cc
               COMPILE_DEFINITIONS ENABLE_CACHING=true
               LABELS unit discretization)
