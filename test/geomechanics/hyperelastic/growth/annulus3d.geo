//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 1.0, 1.0};
//+
Disk(2) = {0, 0, 0, 0.1, 0.1};

BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Extrude {0, 0, 0.1} {
  Surface{1}; Layers {5};
}
