SetFactory("OpenCASCADE");

dx = 500;
Box(1) = {-dx/2.0, -dx/2.0, -dx/2.0, dx, dx, dx};

r = 5;
Disk(7) = {0, 0, 0, r, r};

angle = 20;
Rotate {{0, 1, 0}, {0, 0, 0}, angle*Pi/180.0} { Surface{7}; }

Surface{7} In Volume{1};
Physical Volume(1) = {1};
Physical Surface(1) = {7};

numElementsDomain = 5;
numElementsFracture = 10;
Characteristic Length{1:8} = dx/numElementsDomain;
Characteristic Length{9  } = r*2.0/numElementsFracture;
