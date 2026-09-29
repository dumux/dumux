# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_shallowwater_longwave_gauges
               SOURCES test_longwave_gauges.cc
               LABELS shallowwater unit
               CMAKE_GUARD dune-foamgrid_FOUND)
