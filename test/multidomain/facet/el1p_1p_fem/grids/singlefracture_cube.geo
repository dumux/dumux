a = 1;             // width of the domain
b = 1;             // height of the domain

numCells = 75;
dx_m = (a+b)/2.0/numCells; // discretization length in the matrix
dx_f = dx_m*0.5;       // discretization length in the fracture

h = DefineNumber[1/numCells, Min 0.001, Max 0.1, Step 0.001, Name "Discretization length parameter"];

// Points
Point(1) = {0.0,   0.0, 0.0, dx_m};
Point(2) = {a,     0.0, 0.0, dx_m};
Point(3) = {a,     b,   0.0, dx_m};
Point(4) = {0.0,   b,   0.0, dx_m};
Point(5) = {0.5,   0.0, 0.0, dx_m};
Point(6) = {0.5,   b,   0.0, dx_m};

Point(7) = {a/2.0, 0.25*b, 0.0, dx_f};
Point(8) = {a/2.0, 0.75*b, 0.0, dx_f};

Point(9) = {0.0, 0.25*b, 0.0, dx_m};
Point(10) = {0.0, 0.75*b, 0.0, dx_m};

Point(11) = {a, 0.25*b, 0.0, dx_m};
Point(12) = {a, 0.75*b, 0.0, dx_m};

// domain outline
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 11};
Line(4) = {11, 12};
Line(5) = {12, 3};
Line(6) = {3, 6};
Line(7) = {6, 4};
Line(8) = {4, 10};
Line(9) = {10, 9};
Line(10) = {9, 1};



Line(11) = {8, 7};
//+
Line(12) = {8, 6};
//+
Line(13) = {7, 5};
//+
Line(14) = {9, 7};
//+
Line(15) = {7, 11};
//+
Line(16) = {12, 8};
//+
Line(17) = {8, 10};
//+
Curve Loop(1) = {1, -13, -14, 10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14, -11, 17, 9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {17, -8, -7, -12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {16, 12, -6, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 16, 11, 15};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {15, -3, -2, -13};
//+
Plane Surface(6) = {6};

Transfinite Line{7, 17, 14, 1, 6, 16, 15, 2} = 0.5/h + 1;
Transfinite Line{8, 10, 12, 13, 3, 5} = 0.25/h + 1;
Transfinite Line{11, 4, 9} = 0.5/h + 1;
Physical Line(1) = {11};//, 12, 13};
Physical Surface(1) = {1:6};

Transfinite Surface "*";
Recombine Surface "*";
