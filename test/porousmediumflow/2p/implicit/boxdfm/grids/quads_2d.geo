lc = 0.1;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {20.0, 0.0, 0.0, lc};
Point(3) = {20.0, 20.0, 0.0, lc};
Point(4) = {0.0, 20.0, 0.0, lc};

Point(5) = {0.0, 10.0, 0.0, lc};
Point(6) = {20.0, 10.0, 0.0, lc};
Point(7) = {10.0, 10.0, 0.0, lc};
Point(8) = {10.0, 0.0, 0.0, lc};
Point(9) = {10.0, 20.0, 0.0, lc};

Line(1) = {4, 5};
Line(2) = {5, 1};
Line(3) = {1, 8};
Line(4) = {8, 2};
Line(5) = {2, 6};
Line(6) = {6, 3};
Line(7) = {3, 9};
Line(8) = {9, 4};
Line(9) = {5, 7};
Line(10) = {7, 6};
Line(11) = {9, 7};
Line(12) = {7, 8};

// mesh fracture line
Physical Line(1) = {9};
Transfinite Line{1:12} = 10;

Line Loop(1) = {1, 9, -11, 8};
Plane Surface(1) = {1};
Line Loop(2) = {11, 10, 6, 7};
Plane Surface(2) = {2};
Line Loop(3) = {10, -5, -4, -12};
Plane Surface(3) = {3};
Line Loop(4) = {12, -3, -2, 9};
Plane Surface(4) = {4};

Physical Surface(1) = {1, 2, 3, 4};

// use quads
Transfinite Surface "*";
Recombine Surface "*";
