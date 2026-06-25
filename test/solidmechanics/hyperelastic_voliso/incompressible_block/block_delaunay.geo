SetFactory("OpenCASCADE");

// ------------------------------------------------------------
// Mesh size
// Override from command line with:
// gmsh -setnumber lc 6.25 mesh.geo -3 -format msh2 -o mesh.msh
// ------------------------------------------------------------
DefineConstant[
  lc = {12.5, Min 0.001, Max 50, Step 0.1, Name "lc"}
];

Printf("Using lc = %g", lc);

// ------------------------------------------------------------
// Geometry
// ------------------------------------------------------------
L = 50.0;

// Outer cube: [0,50] x [0,50] x [0,50]
Box(1) = {0, 0, 0, L, L, L};

// Internal box: [25,50] x [25,50] x [0,50]
Box(2) = {25, 25, 0, 25, 25, 50};

// Make internal box faces conforming mesh interfaces
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// ------------------------------------------------------------
// Mesh settings
// ------------------------------------------------------------
Mesh.Algorithm3D = 1; // Delaunay

Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

// Set mesh size on all geometry points
Characteristic Length{ PointsOf{ Volume{:}; } } = lc;

// ------------------------------------------------------------
// Physical groups
// ------------------------------------------------------------
Physical Volume("domain") = Volume{:};
Physical Surface("surfaces") = Surface{:};
