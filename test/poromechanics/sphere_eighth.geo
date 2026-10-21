// One eighth of the sphere of radius R0 for Cryer's benchmark
// Physical surfaces: 1, 2, 3 = symmetry planes x = 0, y = 0, z = 0; 4 = sphere surface
// gmsh -3 -format msh2 -setnumber lc <maximum element size> sphere_eighth.geo -o <mesh file>

SetFactory("OpenCASCADE");

R0 = 1.0;
DefineConstant[ lc = 0.095 ];

Sphere(1) = {0, 0, 0, R0};
Box(2) = {0, 0, 0, 2*R0, 2*R0, 2*R0};
BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; }

Physical Surface("symmetry_z", 3) = {1};
Physical Surface("symmetry_y", 2) = {2};
Physical Surface("symmetry_x", 1) = {3};
Physical Surface("sphere_surface", 4) = {4};
Physical Volume("domain", 1) = {1};

Mesh.MeshSizeMax = lc;
Mesh.Algorithm = 5;
Mesh.Optimize = 1;
