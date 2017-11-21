size = 0.1;
Point(1) = {0.0, 0.0, 0.0, size};
Point(2) = {0.0, 0.0, 0.5, size};

Point(3) = {0.0, 0.0, 1.0, size};
Point(4) = {1.0, 0.0, 1.0, size};
Point(5) = {2.0, 0.0, 1.0, size};
Point(6) = {3.0, -1.0, 1.0, size};

Point(7) = {0.0, 0.0, 0.0, size};
Point(8) = {1.0, 0.0, 0.0, size};
Point(9) = {2.0, 0.0, 0.0, size};
Point(10) = {3.0, -1.0, 0.0, size};

Point(11) = {0.0, -1.0, 1.0, size};
Point(12) = {2.0, 1.0, 1.0, size};
Point(13) = {0.0, -1.0, 0.0, size};
Point(14) = {2.0, 1.0, 0.0, size};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(5) = {11, 13};
Line(8) = {1, 8};
Line(9) = {8, 4};
Line(10) = {4, 3};
Line(11) = {13, 8};
Line(12) = {8, 14};
Line(13) = {14, 12};
Line(14) = {12, 4};
Line(15) = {11, 4};
Line(16) = {4, 5};
Line(17) = {5, 6};
Line(18) = {6, 10};
Line(19) = {10, 9};
Line(20) = {9, 8};
Line(28) = {5, 9};

Line Loop(21) = {10, -2, -1, 8, 9};
Plane Surface(22) = {21};
Line Loop(23) = {14, -9, 12, 13};
Plane Surface(24) = {23};
Line Loop(25) = {15, -9, -11, -5};
Plane Surface(26) = {25};
Line Loop(29) = {17, 18, 19, -28};
Plane Surface(30) = {29};
Line Loop(31) = {16, 28, 20, 9};
Plane Surface(32) = {31};
