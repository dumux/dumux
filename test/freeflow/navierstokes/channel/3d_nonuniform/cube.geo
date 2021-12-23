//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1e-2, 0.5e-2, 0.5e-2};
//+
Physical Surface(1) = {1};
//+
Physical Surface(2) = {2};
//+
Physical Surface(0) = {6, 4, 5, 3};
//+
Physical Volume(0) = {1};

Transfinite Line "*" = 3;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
