# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "delaunay3d.msh")

dumux_add_test(NAME test_boxfvgeometry
              SOURCES test_boxfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false
              LABELS unit discretization)

dumux_add_test(NAME test_boxfvgeometry_caching
              SOURCES test_boxfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true
              LABELS unit discretization)

dumux_add_test(NAME test_boxvolumes
              SOURCES test_boxvolumes.cc
              LABELS unit discretization)