# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_2pve_columnmapping
               SOURCES test_columnmapping.cc
               LABELS unit porousmediumflow 2pve)

dumux_add_test(NAME test_2pve_spatialparams
               SOURCES test_spatialparams.cc
               LABELS unit porousmediumflow 2pve)

dumux_add_test(NAME test_2pve_finelevelview
               SOURCES test_finelevelview.cc
               LABELS unit porousmediumflow 2pve)

dumux_add_test(NAME test_2pve_quantityreconstruction
               SOURCES test_quantityreconstruction.cc
               LABELS unit porousmediumflow 2pve)
