lc = 7.5;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {50.0, 0.0, 0.0, lc};
Point(3) = {50.0, 50.0, 0.0, lc};
Point(4) = {0.0, 50.0, 0.0, lc};

Point(5) = {0.0, 0.0, 25.0, lc};
Point(6) = {40.0, 0.0, 25.0, lc};
Point(7) = {40.0, 50.0, 25.0, lc};
Point(8) = {0.0, 50.0, 25.0, lc};

Point(9) = {0.0, 0.0, 50.0, lc};
Point(10) = {50.0, 0.0, 50.0, lc};
Point(11) = {50.0, 50.0, 50.0, lc};
Point(12) = {0.0, 50.0, 50.0, lc};

Line(1) = {5, 6};
Line(2) = {6, 7};
Line(3) = {7, 8};
Line(4) = {8, 5};
Line(5) = {5, 1};
Line(6) = {1, 4};
Line(7) = {4, 8};
Line(8) = {4, 3};
Line(9) = {3, 2};
Line(10) = {2, 1};
Line(11) = {2, 6};
Line(12) = {7, 3};
Line(13) = {3, 11};
Line(14) = {11, 10};
Line(15) = {10, 2};
Line(16) = {10, 6};
Line(17) = {7, 11};
Line(18) = {11, 12};
Line(19) = {12, 9};
Line(20) = {9, 5};
Line(21) = {8, 12};
Line(22) = {9, 10};

// explicitly mesh fracture surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Line Loop(2) = {4, -20, -19, -21};
Plane Surface(2) = {2};
Line Loop(3) = {21, -18, -17, 3};
Plane Surface(3) = {3};
Line Loop(4) = {16, -1, -20, 22};
Plane Surface(4) = {4};
Line Loop(5) = {17, 14, 16, 2};
Plane Surface(5) = {5};
Line Loop(6) = {14, -22, -19, -18};
Plane Surface(6) = {6};
Surface Loop(1) = {4, 5, 3, 2, 1, 6};
Volume(1) = {1};
Line Loop(7) = {13, -17, 12};
Plane Surface(7) = {7};
Line Loop(8) = {12, 9, 11, 2};
Plane Surface(8) = {8};
Line Loop(9) = {9, -15, -14, -13};
Plane Surface(9) = {9};
Line Loop(10) = {15, 11, -16};
Plane Surface(10) = {10};
Line Loop(11) = {11, -1, 5, -10};
Plane Surface(11) = {11};
Line Loop(12) = {10, 6, 8, 9};
Plane Surface(12) = {12};
Line Loop(13) = {6, 7, 4, 5};
Plane Surface(13) = {13};
Line Loop(14) = {7, -3, 12, -8};
Plane Surface(14) = {14};
Surface Loop(2) = {10, 9, 7, 5, 8};
Volume(2) = {2};
Surface Loop(3) = {8, 1, 13, 12, 11, 14};
Volume(3) = {3};
Physical Volume(1) = {1, 2, 3};
