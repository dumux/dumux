# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_md_mandel_analytical_solution
              LABELS multidomain poromechanics 1p poroelastic
              SOURCES test_analyticalsolution.cc
              CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../params.input)
