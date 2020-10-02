SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
Sphere(2) = {0, 0, 0, 0.7, -Pi/2, Pi/2, 2*Pi};
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

Mesh.MshFileVersion = 2.0;
Mesh.CharacteristicLengthMax = 0.1;
