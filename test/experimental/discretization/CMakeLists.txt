# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_disc_fvgridvariables
               LABELS unit experimental
               SOURCES test_fvgridvariables.cc
               COMMAND ./test_disc_fvgridvariables
               CMD_ARGS -Problem.Name gridvarstest)
