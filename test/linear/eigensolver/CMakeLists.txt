# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_eigensolver.cc
               LABELS linear unit
               CMAKE_GUARD "( DUMUX_HAVE_SPECTRA AND HAVE_UMFPACK )")
