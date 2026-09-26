SetFactory("OpenCASCADE");

DefineConstant[ L = 1.0 ];
Rectangle(1) = {0, 0, 0, L, L};
Rectangle(2) = {0.5*L, 0.5*L, 0, 0.5*L, 0.5*L};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
