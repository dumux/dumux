size = 0.7;
length = 10.0;
radius = 1.0;
a = 30/360*2*Pi;

// quadratic elements for boundary parametrization
Mesh.ElementOrder = 2;

Point(1) = {0.0, 0.0, 0.0, size};
Point(17) = {0.0, 0.0, radius, size};
Point(18) = {0.0, 0.0, -radius, size};

Point(2) = {0.0, length, 0.0, size};
Point(3) = {radius, length, 0.0, size};
Point(4) = {-radius, length, 0.0, size};
Point(5) = {0.0, length, radius, size};
Point(6) = {0.0, length, -radius, size};

Point(7) = {length*Cos(a), -length*Sin(a), 0.0, size};
Point(8) = {length*Cos(a) + radius*Sin(a), -length*Sin(a) + radius*Cos(a), 0.0, size};
Point(9) = {length*Cos(a) - radius*Sin(a), -length*Sin(a) - radius*Cos(a), 0.0, size};
Point(10) = {length*Cos(a), -length*Sin(a), radius, size};
Point(11) = {length*Cos(a), -length*Sin(a), -radius, size};

Point(12) = {-length*Cos(a), -length*Sin(a), 0.0, size};
Point(13) = {-length*Cos(a) - radius*Sin(a), -length*Sin(a) + radius*Cos(a), 0.0, size};
Point(14) = {-length*Cos(a) + radius*Sin(a), -length*Sin(a) - radius*Cos(a), 0.0, size};
Point(15) = {-length*Cos(a), -length*Sin(a), radius, size};
Point(16) = {-length*Cos(a), -length*Sin(a), -radius, size};

Point(19) = {radius, 2.0*radius*Cos(a)/3.0, 0.0, size};
Point(20) = {-radius, 2.0*radius*Cos(a)/3.0, 0.0, size};
Point(21) = {0.0, -4.0*radius*Cos(a)/3.0, 0.0, size};

//optional circle segments
//Point(22) = {0, 2.0*radius*Cos(a)/3.0, 0.0, size};
//Point(23) = {0, 2.0*radius*Cos(a)/3.0, radius, size};
//Point(24) = {0, 2.0*radius*Cos(a)/3.0, -radius, size};

Circle(1) = {15, 12, 14};
Circle(2) = {14, 12, 16};
Circle(3) = {16, 12, 13};
Circle(4) = {13, 12, 15};
Circle(5) = {4, 2, 6};
Circle(6) = {4, 2, 5};
Circle(7) = {5, 2, 3};
Circle(8) = {6, 2, 3};
Circle(9) = {8, 7, 10};
Circle(10) = {10, 7, 9};
Circle(11) = {9, 7, 11};
Circle(12) = {11, 7, 8};
Line(13) = {15, 17};
Line(15) = {17, 10};
Line(17) = {17, 5};
Line(18) = {6, 18};
Line(19) = {18, 11};
Line(20) = {16, 18};
Line(21) = {13, 20};
Line(22) = {4, 20};
Line(23) = {3, 19};
Line(24) = {19, 8};
Line(25) = {21, 9};
Line(26) = {21, 14};
Ellipse(27) = {17, 1, 18, 21};
Ellipse(28) = {18, 1, 17, 21};
Ellipse(29) = {17, 1, 18, 19};
Ellipse(30) = {18, 1, 17, 19};
Ellipse(31) = {17, 1, 18, 20};
Ellipse(32) = {18, 1, 17, 20};

//Circle(33) = {24, 22, 19};
//Circle(34) = {19, 22, 23};
//Circle(35) = {23, 22, 20};
//Circle(36) = {20, 22, 24};
Line Loop(33) = {5, 8, -7, -6};
Plane Surface(34) = {33};
Line Loop(35) = {11, 12, 9, 10};
Plane Surface(36) = {35};
Line Loop(37) = {3, 4, 1, 2};
Plane Surface(38) = {37};
Line Loop(39) = {18, 32, -22, 5};
Ruled Surface(40) = {39};
Line Loop(41) = {32, -21, -3, 20};
Ruled Surface(42) = {41};
Line Loop(43) = {18, 30, -23, -8};
Ruled Surface(44) = {43};
Line Loop(45) = {22, -31, 17, -6};
Ruled Surface(46) = {45};
Line Loop(47) = {7, 23, -29, 17};
Ruled Surface(48) = {47};
Line Loop(49) = {21, -31, -13, -4};
Ruled Surface(50) = {49};
Line Loop(51) = {13, 27, 26, -1};
Ruled Surface(52) = {51};
Line Loop(53) = {20, 28, 26, 2};
Ruled Surface(54) = {53};
Line Loop(55) = {25, -10, -15, 27};
Ruled Surface(56) = {55};
Line Loop(57) = {25, 11, -19, 28};
Ruled Surface(58) = {57};
Line Loop(59) = {15, -9, -24, -29};
Ruled Surface(60) = {59};
Line Loop(61) = {19, 12, -24, -30};
Ruled Surface(62) = {61};

Line(64) = {17, 1};
Line(65) = {1, 18};
Line Loop(66) = {29, -30, -65, -64};
Plane Surface(67) = {66};
Line Loop(68) = {28, -27, 64, 65};
Plane Surface(69) = {68};

Line Loop(70) = {31, -32, -65, -64};
Plane Surface(71) = {70};
Surface Loop(72) = {34, 40, 44, 48, 46, 67, 71};
Volume(73) = {72};
Surface Loop(74) = {42, 50, 52, 54, 38, 71, 69};
Volume(75) = {74};
Surface Loop(76) = {58, 56, 36, 62, 60, 69, 67};
Volume(77) = {76};

Physical Surface(1) = {34};
Physical Surface(2) = {36};
Physical Surface(3) = {38};
Physical Surface(4) = {40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62};
Physical Volume(1) = {73};
Physical Volume(2) = {75};
Physical Volume(3) = {77};
