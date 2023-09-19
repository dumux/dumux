//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 1.0, 1.0};
//+
Disk(2) = {0, 0, 0, 0.1, 0.1};
//+
Disk(3) = {0, 0, 0, 0.95, 0.95};

BooleanDifference{ Surface{3}; Delete; }{ Surface{2}; }
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; }

Field[1] = MathEval;
Field[1].F = "max(0.5 - 0.5*(sqrt(x*x + y*y)+0.05), 0.01)";

Background Field = 1;
