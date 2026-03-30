SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0.0, 0, 1, 1, 0};
Point(5) = {0.25, 0.25, 0, 0.05};
Point(6) = {0.75, 0.75, 0, 0.05};
Line(5) = {5, 6};
Line{5} In Surface{1};
Physical Line(1) = {5};
Physical Surface(1) = {1};
