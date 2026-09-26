# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_diffusioncontainers.cc
               LABELS unit flux)

dumux_add_test(NAME test_tracefields
               SOURCES test_tracefields.cc
               LABELS unit flux)
