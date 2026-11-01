# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_1p_pointsources_interface
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS ADDPOINTSOURCES=0)

dumux_add_test(NAME test_1p_pointsources_interface_addpointsources
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS ADDPOINTSOURCES=1
              EXPECT_COMPILE_FAIL)

dumux_add_test(NAME test_1p_pointsources_interface_addpointsources_template
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS ADDPOINTSOURCES=2
              EXPECT_COMPILE_FAIL)
