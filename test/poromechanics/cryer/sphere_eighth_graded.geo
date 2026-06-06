// 1/8 sphere mesh for the Cryer consolidation benchmark — GRADED variant.
// Radius R0 = 1 m, symmetry planes at x=0, y=0, z=0.
//
// The Cryer/Mandel-Cryer physics has the sharpest gradients at the drained
// outer surface (the consolidation/drainage front), while the field at the
// centre (where the benchmark probes the pressure) is smooth. We therefore
// grade the mesh: fine at the sphere surface, coarse toward the centre.
//
// Physical boundary tags (used in DuMux BCs via GridManager):
//   1 = x=0 symmetry plane,  2 = y=0,  3 = z=0,  4 = sphere surface (drained)
//
// Tune via command line, e.g.:
//   gmsh -3 -format msh2 sphere_eighth_graded.geo -o sphere_eighth.msh \
//        -setnumber lcMin 0.04 -setnumber lcMax 0.15 -setnumber distMax 0.5

SetFactory("OpenCASCADE");

R0 = 1.0;
If (!Exists(lcMin))   lcMin   = 0.04; EndIf  // target size at the sphere surface
If (!Exists(lcMax))   lcMax   = 0.15; EndIf  // target size in the interior/centre
If (!Exists(distMax)) distMax = 0.5;  EndIf  // grading transition distance (in units of R0)

// Build 1/8 sphere by intersecting the full sphere with the positive-octant box
Sphere(1) = {0, 0, 0, R0};
Box(2)    = {0, 0, 0, 2*R0, 2*R0, 2*R0};
BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// OpenCASCADE ordering for Sphere ∩ positive-octant box:
//   Surface 1 = z=0 disk, 2 = y=0 disk, 3 = x=0 disk, 4 = sphere quarter
Physical Surface("symmetry_z", 3) = {1};
Physical Surface("symmetry_y", 2) = {2};
Physical Surface("symmetry_x", 1) = {3};
Physical Surface("sphere_surface", 4) = {4};
Physical Volume("domain", 1) = {1};

// Graded background field: size = lcMin on the sphere surface, growing to lcMax
// over a distance distMax away from it.
Field[1] = Distance;
Field[1].SurfacesList = {4};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcMin;
Field[2].SizeMax = lcMax;
Field[2].DistMin = 0.0;
Field[2].DistMax = distMax;
Background Field = 2;

// Let the background field fully control sizing
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.Algorithm   = 5;  // Delaunay 2D
Mesh.Algorithm3D = 4;  // Frontal 3D
Mesh.Optimize    = 1;
