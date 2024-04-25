// domain measures
xmax = 1.0;
ymax = 1.0;
mesh_size = 0.05;

Point(0) = {0, 0, 0, mesh_size};
Point(1) = {xmax, 0, 0, mesh_size};
Point(2) = {xmax, ymax, 0, mesh_size};
Point(3) = {0, ymax, 0, mesh_size};

Line(0) = {0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 0};

Curve Loop(1) = {0, 1, 2, 3};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
Physical Curve(1) = {1};
