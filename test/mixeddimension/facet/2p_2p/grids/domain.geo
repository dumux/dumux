size = 30;
height = 150;

// fracture corners and intersection points (bottom layer)
Point(1) = {44.7964, 528.5974, 0, size};
Point(2) = {92.08839999999999, 501.1482, 0, size};
Point(3) = {270.148, 397.7989, 0, size};
Point(4) = {356.0905, 347.6924, 0, size};
Point(5) = {270.6623, 524.7026, 0, size};
Point(6) = {269.1347, 129.9814, 0, size};
Point(7) = {32.2042, 309.3757, 0, size};
Point(8) = {25.1027, 285.1246, 0, size};
Point(9) = {0, 333.7332, 0, size};
Point(10) = {441.2444, 0, 0, size};
Point(11) = {364.5808, 24.3622, 0, size};
Point(12) = {387.6907, 40.5051, 0, size};
Point(13) = {700, 275.1918, 0, size};
Point(14) = {269.1347, 79.9814, 0, size};
Point(15) = {107.0347, 543.9814, 0, size};

// fracture corners and intersection points (top layer)
Point(16) = {44.7964, 528.5974, height, size};
Point(17) = {92.08839999999999, 501.1482, height, size};
Point(18) = {270.148, 397.7989, height, size};
Point(19) = {356.0905, 347.6924, height, size};
Point(20) = {270.6623, 524.7026, height, size};
Point(21) = {269.1347, 129.9814, height, size};
Point(22) = {32.2042, 309.3757, height, size};
Point(23) = {25.1027, 285.1246, height, size};
Point(24) = {0, 333.7332, height, size};
Point(25) = {441.2444, 0, height, size};
Point(26) = {364.5808, 24.3622, height, size};
Point(27) = {387.6907, 40.5051, height, size};
Point(28) = {700, 275.1918, height, size};
Point(29) = {269.1347, 79.9814, height, size};
Point(30) = {107.0347, 543.9814, height, size};

// domain corners
Point(31) = {0, 0, 0, size};
Point(32) = {700, 0, 0, size};
Point(33) = {700, 600, 0, size};
Point(34) = {0, 600, 0, size};
Point(35) = {0, 0, height, size};
Point(36) = {700, 0, height, size};
Point(37) = {700, 600, height, size};
Point(38) = {0, 600, height, size};

// domain outline
Line(1) = {38, 24};
Line(2) = {24, 35};
Line(3) = {35, 25};
Line(4) = {25, 36};
Line(5) = {36, 28};
Line(6) = {28, 37};
Line(7) = {37, 38};
Line(8) = {38, 34};
Line(9) = {37, 33};
Line(10) = {36, 32};
Line(11) = {35, 31};
Line(12) = {31, 10};
Line(13) = {10, 32};
Line(14) = {32, 13};
Line(15) = {13, 33};
Line(16) = {33, 34};
Line(17) = {34, 9};
Line(18) = {9, 31};

// domain faces
Line Loop(19) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(20) = {19};
Line Loop(21) = {5, 6, 9, -15, -14, -10};
Plane Surface(22) = {21};
Line Loop(23) = {16, -8, -7, 9};
Plane Surface(24) = {23};
Line Loop(25) = {8, 17, 18, -11, -2, -1};
Plane Surface(26) = {25};
Line Loop(27) = {16, 17, 18, 12, 13, 14, 15};
Plane Surface(28) = {27};
Line Loop(29) = {10, -13, -12, -11, 3, 4};
Plane Surface(30) = {29};

// the fracture edges on the top face
Line(31) = {24, 22};
Line(32) = {23, 22};
Line(33) = {22, 17};
Line(34) = {17, 30};
Line(35) = {16, 17};
Line(36) = {17, 18};
Line(37) = {18, 19};
Line(38) = {20, 18};
Line(39) = {18, 21};
Line(40) = {21, 29};
Line(41) = {22, 21};
Line(42) = {21, 27};
Line(43) = {27, 25};
Line(44) = {26, 27};
Line(45) = {27, 28};
Physical Line(1) = {31:45};
Line{31:45} In Surface{20};

// the fracture edge on the left face
Line(46) = {24, 9};
Physical Line(2) = {46};
Line{46} In Surface{26};

// the fracture edge on the front face
Line(47) = {25, 10};
Physical Line(3) = {47};
Line{47} In Surface{30};

// the fracture edge on the right face
Line(48) = {28, 13};
Physical Line(4) = {48};
Line{48} In Surface{22};

// the fracture edges on the lower face
Line(49) = {9, 7};
Line(50) = {7, 8};
Line(51) = {7, 2};
Line(52) = {2, 15};
Line(53) = {1, 2};
Line(54) = {2, 3};
Line(55) = {3, 4};
Line(56) = {3, 5};
Line(57) = {3, 6};
Line(58) = {6, 14};
Line(59) = {7, 6};
Line(60) = {6, 12};
Line(61) = {12, 10};
Line(62) = {11, 12};
Line(63) = {12, 13};
Physical Line(5) = {49:63};
Line{49:63} In Surface{28};

// remaining fracture edges
Line(64) = {23, 8};
Line(65) = {17, 2};
Line(66) = {30, 15};
Line(67) = {16, 1};
Line(68) = {20, 5};
Line(69) = {18, 3};
Line(70) = {19, 4};
Line(71) = {22, 7};
Line(72) = {21, 6};
Line(73) = {29, 14};
Line(74) = {26, 11};
Line(75) = {27, 12};

// fracture faces
Line Loop(76) = {71, 50, -64, 32};
Plane Surface(77) = {76};
Line Loop(78) = {46, 49, -71, -31};
Plane Surface(79) = {78};
Line Loop(80) = {33, 65, -51, -71};
Plane Surface(81) = {80};
Line Loop(82) = {53, -65, -35, 67};
Plane Surface(83) = {82};
Line Loop(84) = {65, 52, -66, -34};
Plane Surface(85) = {84};
Line Loop(86) = {65, 54, -69, -36};
Plane Surface(87) = {86};
Line Loop(88) = {68, -56, -69, -38};
Plane Surface(89) = {88};
Line Loop(90) = {70, -55, -69, 37};
Plane Surface(91) = {90};
Line Loop(92) = {69, 57, -72, -39};
Plane Surface(93) = {92};
Line Loop(94) = {40, 73, -58, -72};
Plane Surface(95) = {94};
Line Loop(96) = {42, 75, -60, -72};
Plane Surface(97) = {96};
Line Loop(98) = {44, 75, -62, -74};
Plane Surface(99) = {98};
Line Loop(100) = {47, -61, -75, 43};
Plane Surface(101) = {100};
Line Loop(102) = {45, 48, -63, -75};
Plane Surface(103) = {102};
Line Loop(104) = {41, 72, -59, -71};
Plane Surface(105) = {104};

Physical Surface(1) = {77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105};

// the domain volume
Surface Loop(104) = {20, 26, 24, 28, 30, 22};
Volume(105) = {104};
Physical Volume(1) = {105};

Surface{77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105} In Volume{105};
