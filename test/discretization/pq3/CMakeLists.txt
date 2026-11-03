# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_pq3dofconsistency
               SOURCES test_pq3dofconsistency.cc
               CMAKE_GUARD dune-alugrid_FOUND
               LABELS unit discretization)
