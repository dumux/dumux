# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_function_l2norm.cc
              CMAKE_GUARD dune-functions_FOUND
              LABELS unit)
