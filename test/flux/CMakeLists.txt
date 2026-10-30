# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(cvfe)

dumux_add_test(SOURCES test_diffusioncontainers.cc
               LABELS unit flux)
