nPoints = 15;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {0.0, 0.25, 0.25};
Point(3) = {0.0, 0.75, 0.75};
Point(4) = {0.0, 1.0,  1.0};

Point(5) = {0.5, 0.0, 0.0};
Point(6) = {0.5, 0.25, 0.25};
Point(7) = {0.5, 0.75, 0.75};
Point(8) = {0.5, 1.0,  1.0};

Point(9) = {1.0, 0.0, 0.0};
Point(10) = {1.0, 0.25, 0.25};
Point(11) = {1.0, 0.75, 0.75};
Point(12) = {1.0, 1.0,  1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Line(4) = {1, 5};
Line(5) = {2, 6};
Line(6) = {3, 7};
Line(7) = {4, 8};

Line(8) = {5, 6};
Line(9) = {6, 7};
Line(10) = {7, 8};

Line(11) = {5, 9};
Line(12) = {6, 10};
Line(13) = {7, 11};
Line(14) = {8, 12};

Line(15) = {9, 10};
Line(16) = {10, 11};
Line(17) = {11, 12};

Curve Loop(1) = {4, 8, -5, -1};
Plane Surface(1) = {1};
Curve Loop(2) = {11, 15, -12, -8};
Plane Surface(2) = {2};
Curve Loop(3) = {12, 16, -13, -9};
Plane Surface(3) = {3};
Curve Loop(4) = {5, 9, -6, -2};
Plane Surface(4) = {4};
Curve Loop(5) = {13, 17, -14, -10};
Plane Surface(5) = {5};
Curve Loop(6) = {6, 10, -7, -3};
Plane Surface(6) = {6};

Transfinite Line{1,3,8,10,15,17} = nPoints/2.0;
Transfinite Line{2,4:7,9,11:14,16} = nPoints;
Transfinite Surface "*";
Recombine Surface "*";

Physical Line(1) = {9};
Physical Surface(1) = {1,2,3,4,5,6};
