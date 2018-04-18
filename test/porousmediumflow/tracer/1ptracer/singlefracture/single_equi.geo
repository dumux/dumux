///////////////////////////////////////////////////////////////////
// Gmsh geometry file for the domain used in Barlag 1998
// The 3d domain consists of two layers and a 2d fault zone
// embedded in the upper one. Meshing is done using tetrahedra.
///////////////////////////////////////////////////////////////////

// fault zone aperture
a = 1.0;

ref = 1.0;                    // refinement in z-direction towards fracture plane (ref < 1.0)
numPointsZ_layer1 = 5;       // no vertices in z in layer one
numPointsZ_layer2_below = 5; // no vertices in z in layer two below the fault
numPointsZ_layer2_above = 5; // no vertices in z in layer two above the fault
numPointsZ_fault = 5;        // no vertices in z inside the fault
numPointsX = 10;             // no vertices in x direction
numPointsY = 10;             // no vertices in y direction

// domain bounding box
Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {100.0, 0.0, 0.0, 1.0};
Point(3) = {100.0, 100.0, 0.0, 1.0};
Point(4) = {0.0, 100.0, 0.0, 1.0};
Point(5) = {0.0, 0.0, 100.0, 1.0};
Point(6) = {100.0, 0.0, 100.0, 1.0};
Point(7) = {100.0, 100.0, 100.0, 1.0};
Point(8) = {0.0, 100.0, 100.0, 1.0};

// Lower layer boundary points
Point(9) = {0.0, 0.0, 10.0, 1.0};
Point(10) = {100.0, 0.0, 10.0, 1.0};
Point(11) = {100.0, 100.0, 10.0, 1.0};
Point(12) = {0.0, 100.0, 10.0, 1.0};

// fault zone boundary points
Point(13) = {0.0, 0.0, 80.0 - a, 1.0};
Point(14) = {100.0, 0.0, 20.0 - a, 1.0};
Point(15) = {100.0, 100.0, 20.0 - a, 1.0};
Point(16) = {0.0, 100.0, 80.0 - a, 1.0};
Point(17) = {0.0, 0.0, 80.0 + a, 1.0};
Point(18) = {100.0, 0.0, 20.0 + a, 1.0};
Point(19) = {100.0, 100.0, 20.0 + a, 1.0};
Point(20) = {0.0, 100.0, 80.0 + a, 1.0};

// layer one vertical discretization
Line(1) = {1, 9};
Line(2) = {4, 12};
Line(3) = {2, 10};
Line(4) = {3, 11};
Transfinite Line{1:4} = numPointsZ_layer1;

// layer two vertical discretization
Line(5) = {13, 9};
Line(6) = {14, 10};
Line(7) = {15, 11};
Line(8) = {16, 12};
Transfinite Line{5:8} = numPointsZ_layer2_below Using Progression ref;
Line(9) = {20, 8};
Line(10) = {19, 7};
Line(11) = {18, 6};
Line(12) = {17, 5};
Transfinite Line{9:12} = numPointsZ_layer2_above Using Progression ref;

// fault zone vertical discretization
Line(13) = {17, 13};
Line(14) = {18, 14};
Line(15) = {19, 15};
Line(16) = {20, 16};
Transfinite Line{13:16} = numPointsZ_fault;

// discretization in x-direction
Line(17) = {1, 2};
Line(18) = {9, 10};
Line(19) = {4, 3};
Line(20) = {12, 11};
Line(21) = {13, 14};
Line(22) = {16, 15};
Line(23) = {17, 18};
Line(24) = {20, 19};
Line(25) = {5, 6};
Line(26) = {8, 7};
Transfinite Line{17:26} = numPointsX;

// discretization in y-direction
Line(27) = {2, 3};
Line(28) = {11, 10};
Line(29) = {14, 15};
Line(30) = {19, 18};
Line(31) = {6, 7};
Line(32) = {8, 5};
Line(33) = {17, 20};
Line(34) = {16, 13};
Line(35) = {9, 12};
Line(36) = {4, 1};
Transfinite Line{27:36} = numPointsY;

/////////////////////////
// Surfaces & volumes
////////////////////////

// lower layer volume
Line Loop(37) = {17, 3, -18, -1};
Plane Surface(38) = {37};
Line Loop(39) = {27, 4, 28, -3};
Plane Surface(40) = {39};
Line Loop(41) = {4, -20, -2, 19};
Plane Surface(42) = {41};
Line Loop(43) = {2, -35, -1, -36};
Plane Surface(44) = {43};
Line Loop(45) = {35, 20, 28, -18};
Plane Surface(46) = {45};
Line Loop(47) = {27, -19, 36, 17};
Plane Surface(48) = {47};
Surface Loop(49) = {46, 44, 42, 40, 48, 38};
Volume(50) = {49};

// layer 2 below fault
Line Loop(51) = {5, 18, -6, -21};
Plane Surface(52) = {51};
Line Loop(53) = {6, -28, -7, -29};
Plane Surface(54) = {53};
Line Loop(55) = {7, -20, -8, 22};
Plane Surface(56) = {55};
Line Loop(57) = {8, -35, -5, -34};
Plane Surface(58) = {57};
Line Loop(59) = {21, 29, -22, 34};
Plane Surface(60) = {59};
Surface Loop(61) = {60, 52, 58, 56, 54, 46};
Volume(62) = {61};

// layer two above fault
Line Loop(63) = {12, 25, -11, -23};
Plane Surface(64) = {63};
Line Loop(65) = {11, 31, -10, 30};
Plane Surface(66) = {65};
Line Loop(67) = {10, -26, -9, 24};
Plane Surface(68) = {67};
Line Loop(69) = {9, 32, -12, 33};
Plane Surface(70) = {69};
Line Loop(71) = {25, 31, -26, 32};
Plane Surface(72) = {71};
Line Loop(73) = {23, -30, -24, -33};
Plane Surface(74) = {73};
Surface Loop(75) = {72, 64, 70, 68, 66, 74};
Volume(76) = {75};

// fault layer
Line Loop(77) = {13, 21, -14, -23};
Plane Surface(78) = {77};
Line Loop(79) = {14, 29, -15, 30};
Plane Surface(80) = {79};
Line Loop(81) = {15, -22, -16, 24};
Plane Surface(82) = {81};
Line Loop(83) = {16, 34, -13, 33};
Plane Surface(84) = {83};
Surface Loop(85) = {78, 84, 82, 80, 74, 60};
Volume(86) = {85};

// use hexahedrons everywhere
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";

// Give the layers physical entity indices
Physical Volume(1) = {50};     // lower layer
Physical Volume(2) = {62, 76}; // upper layer
Physical Volume(3) = {86};     // fault zone
