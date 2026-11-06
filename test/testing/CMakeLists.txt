# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "test_convergencetest.py")

add_executable(test_convergencetest_quadrature EXCLUDE_FROM_ALL test_convergencetest_quadrature.cc)

dumux_add_test(NAME test_convergencetest
               TARGET test_convergencetest_quadrature
               LABELS unit testing
               COMMAND ./test_convergencetest.py
               CMD_ARGS ${DUMUX_CONVERGENCETEST} ${CMAKE_CURRENT_BINARY_DIR}/test_convergencetest_quadrature)
