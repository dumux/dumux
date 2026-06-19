# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_mumps_sequential
               SOURCES main.cc
               CMAKE_GUARD MUMPS_FOUND
               LABELS unit linear)
