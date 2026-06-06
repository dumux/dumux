// 1/8 sphere mesh for the Cryer consolidation benchmark
// Radius R0 = 1 m, symmetry planes at x=0, y=0, z=0
//
// Physical boundary tags (used in DuMux BCs via GridManager):
//   1 = x=0 symmetry plane (quarter-disk in y-z)
//   2 = y=0 symmetry plane (quarter-disk in x-z)
//   3 = z=0 symmetry plane (quarter-disk in x-y)
//   4 = sphere surface (drained outer boundary)

SetFactory("OpenCASCADE");

R0 = 1.0;
// Characteristic length. lc = 0.095 gives a moderate ~854-node uniform mesh — a
// good accuracy/cost balance for this benchmark (the centre pore pressure needs
// adequate *bulk* resolution; boundary grading does not help, see README/notes).
// Decrease for a convergence study, increase for a fast smoke test.
lc = 0.095;

// Build 1/8 sphere by intersecting full sphere with positive-octant box
Sphere(1) = {0, 0, 0, R0};
Box(2)    = {0, 0, 0, 2*R0, 2*R0, 2*R0};
BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// After the Boolean operation the volume is 1, surfaces are auto-numbered 1..4.
// We identify them: three flat quarter-disks + one quarter-sphere.
// OpenCASCADE numbers from the Box intersection produce:
//   Surfaces 1,2,3 = flat faces (x=0, y=0, z=0 in some order)
//   Surface  4     = curved sphere face
// We loop over all surfaces and tag by bounding box.

// The following Macro approach tags by bounding box centroid:
Macro TagSurfaces
  For s In Surface{:}
    bbox() = BoundingBox Surface{s};
    // bbox = {xmin, ymin, zmin, xmax, ymax, zmax}
    cx = (bbox(0)+bbox(3))/2;
    cy = (bbox(1)+bbox(4))/2;
    cz = (bbox(2)+bbox(5))/2;
    // Flat face at x=0: xmin and xmax both ≈ 0
    If (Fabs(bbox(0)) < 1e-6 && Fabs(bbox(3)) < 1e-6)
      Physical Surface("symmetry_x", 1) += {s};
    EndIf
    // Flat face at y=0
    If (Fabs(bbox(1)) < 1e-6 && Fabs(bbox(4)) < 1e-6)
      Physical Surface("symmetry_y", 2) += {s};
    EndIf
    // Flat face at z=0
    If (Fabs(bbox(2)) < 1e-6 && Fabs(bbox(5)) < 1e-6)
      Physical Surface("symmetry_z", 3) += {s};
    EndIf
  EndFor
EndMacro

// The sphere surface is whatever is left after the three flat faces.
// Tag all surfaces first by the macro, then use a catch-all for the sphere.
// For simplicity, just number them explicitly after inspection.
// OpenCASCADE typical ordering for Sphere ∩ Box (positive octant):
//   Surface 1 = z=0 disk
//   Surface 2 = y=0 disk
//   Surface 3 = x=0 disk
//   Surface 4 = sphere quarter
Physical Surface("symmetry_z", 3) = {1};
Physical Surface("symmetry_y", 2) = {2};
Physical Surface("symmetry_x", 1) = {3};
Physical Surface("sphere_surface", 4) = {4};

Physical Volume("domain", 1) = {1};

// Mesh settings.
// NOTE: for this OpenCASCADE geometry the in-file size controls are not honoured
// reliably (the mesh stays near ~220 nodes regardless). Generate the committed
// moderate mesh by passing the size on the command line instead:
//     gmsh -3 -format msh2 -clmax 0.095 sphere_eighth.geo -o sphere_eighth.msh
// (-clmax 0.095 -> ~854 nodes; 0.12 -> ~224 coarse; 0.08 -> ~1311 fine).
MeshSize{ PointsOf{ Volume{1}; } } = lc;
Mesh.Algorithm   = 5;  // Delaunay 2D
Mesh.Algorithm3D = 4;  // Frontal 3D
Mesh.Optimize    = 1;
