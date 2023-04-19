# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_newton.cc
               LABELS unit nonlinear)
dumux_add_test(NAME test_newton_linesearch
               TARGET test_newton
               COMMAND test_newton
               CMD_ARGS "-Newton.UseLineSearch" "true"
               LABELS unit nonlinear)
dumux_add_test(SOURCES test_newton_parallel.cc
               LABELS unit nonlinear parallel
               COMMAND ${MPIEXEC}
               CMD_ARGS -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_newton_parallel
               CMAKE_GUARD HAVE_MPI)
