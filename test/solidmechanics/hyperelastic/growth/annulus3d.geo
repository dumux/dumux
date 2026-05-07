SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 1.0, 1.0};
Disk(2) = {0, 0, 0, 0.1, 0.1};

BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

Extrude {0, 0, 0.1} {
  Surface{1}; Layers{3};
}

Field[1] = MathEval;
Field[1].F = "max(0.5 - 0.5*(sqrt(x*x + y*y)+0.05), 0.02)";

Background Field = 1;
