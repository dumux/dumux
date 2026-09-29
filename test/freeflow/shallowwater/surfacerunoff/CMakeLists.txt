# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_shallowwater_crosssection
               SOURCES test_crosssection.cc
               LABELS shallowwater unit)

dumux_add_test(NAME test_shallowwater_greenampt
               SOURCES test_greenampt.cc
               LABELS shallowwater unit)

dumux_add_test(NAME test_shallowwater_infiltration
               SOURCES test_infiltration.cc
               LABELS shallowwater unit)

dumux_add_test(NAME test_shallowwater_outletcontrol
               SOURCES test_outletcontrol.cc
               LABELS shallowwater unit)
