# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(boundingboxtree)

dumux_add_test(SOURCES test_0d1d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_0d2d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_0d3d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_1d1d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_1d3d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_1d2d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_2d2d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_2d3d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_3d3d_intersection.cc LABELS unit)
dumux_add_test(SOURCES test_distance.cc LABELS unit)
dumux_add_test(SOURCES test_normal.cc LABELS unit)
dumux_add_test(SOURCES test_volume.cc LABELS unit)
dumux_add_test(SOURCES test_triangulate_convex_hull.cc LABELS unit)
dumux_add_test(SOURCES test_intersectingentity_cartesiangrid.cc LABELS unit)
dumux_add_test(SOURCES test_circlepoints.cc LABELS unit)
dumux_add_test(SOURCES test_cylinderintegration.cc LABELS unit)
dumux_add_test(SOURCES test_makegeometry.cc LABELS unit)
dumux_add_test(SOURCES test_refinementquadraturerule.cc LABELS unit)
dumux_add_test(SOURCES test_volumeaverage.cc LABELS unit CMAKE_GUARD dune-subgrid_FOUND)
dune_symlink_to_source_files(FILES ball.msh)
dumux_add_test(SOURCES test_intersectionentityset.cc
               CMAKE_GUARD dune-alugrid_FOUND
               LABELS unit)
