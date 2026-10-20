# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(NAME test_shallowwater_reconstruction
              SOURCES test_shallowwater_reconstruction.cc
              CMAKE_GUARD dune-uggrid_FOUND
              LABELS shallowwater unit)
