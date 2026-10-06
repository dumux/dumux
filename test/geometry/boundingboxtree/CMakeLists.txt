# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# build the tests for the bounding box tree
dumux_add_test(NAME test_bboxtree_dim1
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=1
              LABELS unit)

dumux_add_test(NAME test_bboxtree_dim2
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=2
              LABELS unit)

dumux_add_test(NAME test_bboxtree_dim3
              SOURCES test_bboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=3
              LABELS unit)

dumux_add_test(NAME test_geometry_fracture
              SOURCES test_geometry_fracture.cc
              CMAKE_GUARD dune-foamgrid_FOUND
              LABELS unit)

dumux_add_test(NAME test_3d3d_grid_intersection
              SOURCES test_3d3d_grid_intersection.cc
              LABELS unit)

dumux_add_test(NAME test_distributedbboxtree_dim2
              SOURCES test_distributedbboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=2
              LABELS unit parallel
              MPI_RANKS 1 2 4
              TIMEOUT 1200)

dumux_add_test(NAME test_distributedbboxtree_dim3
              SOURCES test_distributedbboxtree.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=3
              LABELS unit parallel
              MPI_RANKS 1 2 4
              TIMEOUT 1200)

dumux_add_test(NAME test_distributedglue_dim2
              SOURCES test_distributedglue.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=2
              LABELS unit parallel
              MPI_RANKS 1 2 4
              TIMEOUT 1200)

dumux_add_test(NAME test_distributedglue_dim3
              SOURCES test_distributedglue.cc
              COMPILE_DEFINITIONS WORLD_DIMENSION=3
              LABELS unit parallel
              MPI_RANKS 1 2 4
              TIMEOUT 1200)

# symlink the input file in the build directory
dune_symlink_to_source_files(FILES "network1d.msh" "network2d.msh" "fracture.msh")
