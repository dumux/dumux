size = 0.01;
Point(1) = {0.0, 0.0, 0.0, size};
Point(2) = {1.0, 0.0, 0.0, size};
Point(3) = {1.0, 1.0, 0.0, size};
Point(4) = {0.0, 1.0, 0.0, size};
Point(5) = {0.0, 0.5, 0.0, size};
Point(6) = {0.5, 0.5, 0.0, size};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line(6) = {5, 6};

Line Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Line{6} In Surface{1};

Physical Surface(1) = {1};
Physical Line(1) = {6};
