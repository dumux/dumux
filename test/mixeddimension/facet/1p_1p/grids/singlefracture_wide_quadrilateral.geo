numElements = 160;
domainLength = 1;
lc = 1/numElements;

Point(1) = {1.0, -0.5, 0, lc};
Point(2) = {0, 0.5, 0, lc};
Point(3) = {1.0, 0.5, 0, lc};
Point(4) = {-1.0, 0.5, 0, lc};
Point(5) = {0, -0.5, 0, lc};
Point(6) = {-1.0, -0.5, 0, lc};

Line(1) = {1, 3};
Transfinite Line{1} = numElements/2 + 1;
Line(2) = {3, 2};
Transfinite Line{2} = numElements/2 + 1;
Line(3) = {2, 5};
Transfinite Line{3} = numElements/2 + 1;
Line(4) = {5, 1};
Transfinite Line{4} = numElements/2 + 1;
Line(5) = {5, 6};
Transfinite Line{5} = numElements/2 + 1;
Line(6) = {6, 4};
Transfinite Line{6} = numElements/2 + 1;
Line(7) = {4, 2};
Transfinite Line{7} = numElements/2 + 1;

// we want line 2 to be explicitly meshed
Physical Line(1) = {3};

// the two surfaces
Line Loop(8) = {1, 2, 3, 4};
Plane Surface(9) = {8};
Line Loop(10) = {5, 6, 7, 3};
Plane Surface(11) = {10};

Transfinite Surface{9} = {1, 3, 2, 5};
Transfinite Surface{11} = {5, 2, 4, 6};
Recombine Surface{9};
Recombine Surface{11};

Physical Surface(1) = {9, 11};
