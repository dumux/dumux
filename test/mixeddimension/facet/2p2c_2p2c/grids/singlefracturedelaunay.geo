numElements = 30;
lc = 1/numElements;

Point(1) = {0.5, -0.5, 0, lc};
Point(2) = {0, 0, 0, lc};

Point(3) = {0.5, 0.5, 0, lc};
Point(4) = {-0.5, 0.5, 0, lc};
Point(5) = {-0.5, 0, 0, lc};
Point(6) = {-0.5, -0.5, 0, lc};

Point(7) = {0.5, 0.25, 0, lc};
Point(8) = {0.5, -0.25, 0, lc};

Line(1) = {5, 6};
Line(2) = {6, 1};
Line(3) = {1, 8};
Line(4) = {8, 7};
Line(5) = {7, 3};
Line(6) = {3, 4};
Line(7) = {4, 5};

Line Loop(8) = {7, 1, 2, 3, 4, 5, 6};
Plane Surface(9) = {8};

Line(10) = {5, 2};
Line(11) = {2, 7};
// Line(12) = {2, 8};

// we want the fracture lines to be explicitly meshed
Line{10} In Surface{9};
Physical Line(1) = {10};
Line{11} In Surface{9};
Physical Line(2) = {11};
// Line{12} In Surface{9};
// Physical Line(3) = {12};

Physical Surface(1) = {9};
