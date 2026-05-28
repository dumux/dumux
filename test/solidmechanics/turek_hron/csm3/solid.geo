//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0.2, 0.19, 0, 0.4, 0.02, 0};
//+
Disk(2) = {0.2, 0.2, 0, 0.05, 0.05};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
