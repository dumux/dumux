// Rectangular structured cube grid
// compile with gmsh -2 rectangle-domain.geo

Point(1) = {0, 0, 0, 1.0};
Point(2) = {6, 0, 0, 1.0};
Point(3) = {6, 4, 0, 1.0};
Point(4) = {0, 4, 0, 1.0};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(7) = {6};

Transfinite Line{1, 3} = 49; // 48 cells
Transfinite Line{2, 4} = 33; // 33 cells
Transfinite Surface{7}; // structured surface
Recombine Surface{7}; // make quarilateral
