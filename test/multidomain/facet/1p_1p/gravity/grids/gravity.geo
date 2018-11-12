// domain measures
xmax = 1.0;
ymax = 3.0;

// use zmax != 0 to produce a surface grid
zmax = 0.0;

// discretization specifics
numCellsX = 30;
numCellsY = 90;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {xmax, 0.0, 0.0};
Point(3) = {xmax, ymax/2.0, zmax/2.0};
Point(4) = {xmax, ymax, zmax};
Point(5) = {0.0, ymax, zmax};
Point(6) = {0.0, ymax/2.0, zmax/2.0};

Line(1) = {1, 2}; Transfinite Line{1} = numCellsX;
Line(2) = {2, 3}; Transfinite Line{2} = numCellsY/2;
Line(3) = {3, 4}; Transfinite Line{3} = numCellsY/2;
Line(4) = {4, 5}; Transfinite Line{4} = numCellsX;
Line(5) = {5, 6}; Transfinite Line{5} = numCellsY/2;
Line(6) = {6, 1}; Transfinite Line{6} = numCellsY/2;
Line(7) = {3, 6}; Transfinite Line{7} = numCellsX;

Line Loop(1) = {1, 2, 7, 6};
Line Loop(2) = {3, 4, 5, -7};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Physical Line(1) = {7};
Physical Surface(1) = {1, 2};
Transfinite Surface "*";
Recombine Surface "*";
