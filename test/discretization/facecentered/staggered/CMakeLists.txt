# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_facecentered_staggered
               SOURCES test_staggeredfvgeometry.cc
               LABELS unit discretization)

dumux_add_test(NAME test_velocity_reconstruction
               SOURCES test_velocity_reconstruction.cc
               CMAKE_GUARD dune-spgrid_FOUND
               LABELS unit discretization brinkman)
