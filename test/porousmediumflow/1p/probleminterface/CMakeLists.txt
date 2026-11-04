# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_1p_probleminterface
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS INTERFACE=0)

dumux_add_test(NAME test_1p_probleminterface_neumannatpos
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS INTERFACE=1
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_1p_probleminterface_neumann
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS INTERFACE=2
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_1p_probleminterface_neumann_boundaryflux
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS INTERFACE=3)
