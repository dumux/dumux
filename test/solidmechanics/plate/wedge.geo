SetFactory("OpenCASCADE");

DefineConstant[ L = 1.0 ];
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {0, L, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};
