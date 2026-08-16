// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// The benchmark box [-1,1] x [-1,1] x [0,1], meshed with tetrahedra.
//
// Splitting a hexahedral grid into tetrahedra produces slivers, so a mesh made that way
// would test the mesher as much as the coupling. Gmsh's Delaunay refinement gives
// well-shaped elements, so the unstructured test measures the point location and the source
// distribution instead. The element size is set on the command line via -clmax.
SetFactory("OpenCASCADE");
Box(1) = {-1, -1, 0, 2, 2, 1};
Mesh.Algorithm3D = 1;   // Delaunay
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
