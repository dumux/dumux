// domain measures
xmax = 3.0;
ymax = 1.0;

// use zmax != 0 to produce a surface grid
zmax = 0.0;

// discretization specifics
numCellsX = 120;
numCellsY = 40;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {xmax/2.0, 0.0, zmax/2.0};
Point(3) = {xmax, 0.0, zmax};
Point(4) = {xmax, ymax, zmax};
Point(5) = {xmax/2.0,ymax, zmax/2.0};
Point(6) = {0.0, ymax, 0.0};

Line(1) = {1, 2}; Transfinite Line{1} = numCellsX/2;
Line(2) = {2, 3}; Transfinite Line{2} = numCellsX/2;
Line(3) = {3, 4}; Transfinite Line{3} = numCellsY;
Line(4) = {4, 5}; Transfinite Line{4} = numCellsX/2;
Line(5) = {5, 6}; Transfinite Line{5} = numCellsX/2;
Line(6) = {6, 1}; Transfinite Line{6} = numCellsY;
Line(7) = {2, 5}; Transfinite Line{7} = numCellsY;

Line Loop(1) = {1, 7, 5, 6};
Line Loop(2) = {2, 3, 4, -7};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Physical Line(1) = {7};
Physical Surface(1) = {1, 2};
Transfinite Surface "*";
Recombine Surface "*";
