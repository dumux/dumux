size = 0.2;
Point(1) = {0.0, 0.0, 0.0, size};
Point(2) = {1.0, 0.0, 0.0, size};
Point(3) = {1.0, 1.0, 0.0, size};
Point(4) = {0.0, 1.0, 0.0, size};

Point(5) = {0.5, 0.0, 0.0, size};
Point(6) = {0.5, 1.0, 0.0, size};

Line(1) = {3, 6};
Line(6) = {6, 4};
Line(2) = {4, 1};
Line(3) = {2, 5};
Line(4) = {5, 1};
Line(5) = {2, 3};
Line(7) = {6, 5};

Line Loop(8) = {6, 2, -4, -7};
Plane Surface(9) = {8};
Line Loop(10) = {1, 7, -3, 5};
Plane Surface(11) = {10};

Physical Surface(2) = {9};
Physical Surface(1) = {11};

Recombine Surface{9, 11};
