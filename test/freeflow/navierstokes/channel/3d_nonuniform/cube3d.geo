//+
Point(1) = {0, 0, 0, 1};
//+
Point(2) = {1, 0, 0, 1};
//+
Point(3) = {1, 1, 0, 1};
//+
Point(4) = {0, 1, 0, 1};
//+
Point(5) = {0, 1, 1, 1};
//+
Point(6) = {0, 0, 1, 1};
//+
Point(7) = {1, 0, 1, 1};
//+
Point(8) = {1, 1, 1, 1};
//+
Point(9) = {1, 1, 0.4, 1};
//+
Point(10) = {1, 0, 0.4, 1};
//+
Point(11) = {0, 1, 0.7, 1};
//+
Point(12) = {0, 0, 0.7, 1};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {2, 10};
//+
Line(6) = {10, 9};
//+
Line(7) = {9, 3};
//+
Line(8) = {9, 8};
//+
Line(9) = {8, 7};
//+
Line(10) = {7, 10};
//+
Line(11) = {8, 5};
//+
Line(12) = {5, 6};
//+
Line(13) = {6, 7};
//+
Line(14) = {5, 11};
//+
Line(15) = {11, 12};
//+
Line(16) = {12, 6};
//+
Line(17) = {12, 10};
//+
Line(18) = {11, 9};
//+
Line(19) = {12, 1};
//+
Line(20) = {11, 4};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {19, 2, 5, -17};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {19, -1, -20, 15};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {18, 7, 4, -20};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {7, -3, 5, 6};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {18, 8, 11, 14};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {8, 9, 10, 6};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {10, -17, 16, 13};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {15, 16, -12, 14};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {13, -9, 11, 12};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {17, 6, -18, 15};
//+
Plane Surface(11) = {11};
//+
Surface Loop(1) = {3, 2, 1, 5, 4, 11};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {10, 8, 7, 6, 9, 11};
//+
Volume(2) = {2};
//+
Physical Surface(1) = {4, 6};
//+
Physical Surface(2) = {2, 8};
//+
Physical Surface(0) = {5, 1, 3, 9, 10, 7};
//+
Physical Volume(0) = {1, 2};

Transfinite Line "*" = 3;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
