# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_ff_mass_fe_boundaryflux.cc
               LABELS unit freeflow navierstokes)

dumux_add_test(SOURCES test_ff_mass_density_advection.cc
               LABELS freeflow navierstokes
               CMAKE_GUARD HAVE_UMFPACK)
