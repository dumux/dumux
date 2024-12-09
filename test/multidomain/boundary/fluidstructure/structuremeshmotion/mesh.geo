//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0.2, 0.19, 0, 0.4, 0.02, 0};
//+
Disk(2) = {0.2, 0.2, 0, 0.05, 0.05};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Rectangle(3) = {0, 0, 0, 2.5, 0.4, 0};
//+
Disk(4) = {0.2, 0.2, 0, 0.05, 0.05};
//+
BooleanDifference{ Surface{3}; Delete; }{ Surface{4}; Delete; }
//+
BooleanFragments{ Surface{3}; Delete; }{ Surface{1}; Delete; }
//+
//+
Physical Surface(1) = {1};
//+
Physical Surface(2) = {2};
//+
MeshSize {5, 9, 6, 7, 8} = 0.01;
//+
MeshSize {3, 1, 4, 2} = 0.1;
