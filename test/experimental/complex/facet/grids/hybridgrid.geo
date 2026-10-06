numElemsPerSide = 25;
dx = 1/numElemsPerSide;

Point(1) = {0.5, -0.5, 0, dx};
Point(2) = {0.5, 0, 0, dx};
Point(3) = {0.5, 0.5, 0, dx};
Point(4) = {-0.5, 0.5, 0, dx};
Point(5) = {-0.5, 0, 0, dx};
Point(6) = {-0.5, -0.5, 0, dx};

Line(1) = {4, 5};
Transfinite Line{1} = Round(numElemsPerSide/2) + 1;
Line(2) = {5, 2};
Transfinite Line{2} = numElemsPerSide + 1;
Line(3) = {2, 3};
Transfinite Line{3} = Round(numElemsPerSide/2) + 1;
Line(4) = {3, 4};
Transfinite Line{4} = numElemsPerSide + 1;
Line(5) = {5, 6};
Transfinite Line{5} = Round(numElemsPerSide/2) + 1;
Line(6) = {6, 1};
Transfinite Line{6} = numElemsPerSide + 1;
Line(7) = {1, 2};
Transfinite Line{7} = Round(numElemsPerSide/2) + 1;

// we want line 2 to be explicitly meshed
Physical Line(1) = {2};

// the two surfaces
Line Loop(8) = {2, -7, -6, -5};
Plane Surface(9) = {8};
Line Loop(10) = {2, 3, 4, 1};
Plane Surface(11) = {10};

Transfinite Surface{9} = {1, 2, 5, 6};
Transfinite Surface{11} = {2, 3, 4, 5};
Recombine Surface{9};
Recombine Surface{11};

Physical Surface(1) = {9, 11};
