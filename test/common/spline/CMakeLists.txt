# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_spline.cc
              LABELS unit)
dumux_add_test(SOURCES test_cubicspline.cc
              LABELS unit)
dumux_add_test(SOURCES test_monotonecubicspline.cc
              LABELS unit)
