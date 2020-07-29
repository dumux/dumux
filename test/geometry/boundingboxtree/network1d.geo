cl_ = 0.5;

Point(1) = {0.5, 0.5, 0.5, cl_};
Point(2) = {0.5, 0.5, 0.9, cl_};
Point(3) = {0.1, 0.1, 0.1, cl_};
Point(4) = {0.1, 0.9, 0.1, cl_};

Point(5) = {0.9, 0.9, 0.1, cl_};
Point(6) = {0.9, 0.9, 0.9, cl_};

Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {1, 4};
Line(4) = {5, 6};

Transfinite Line{1, 2, 3} = 1;
Transfinite Line{4} = 2;
