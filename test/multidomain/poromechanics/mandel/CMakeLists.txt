# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(analyticalsolution)

dumux_add_test(NAME test_md_mandel
               SOURCES main.cc
               LABELS benchmark mandel 2d)
