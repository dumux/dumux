length = 1.0;
depth = 0.2;
numElementsX = 100;

// compute dummy element size to define at points
size = length/numElementsX;

Point(1) = {0.0, 0.0, 0.0, size};
Point(2) = {length/5.0, -depth/3.0, 0.0, size};
Point(3) = {2*length/5.0, -depth/25.0, 0.0, size};
Point(4) = {4*length/5.0, -depth, 0.0, size};
Point(5) = {length, 0.0, 0.0, size};

// top (boundary) line
Line(1) = {5, 1};

// bottom line
Spline(2) = {5, 4, 3, 2, 1};

// define domain surface
Curve Loop(1) = {2, -1};
Plane Surface(1) = {1};

// fix number of elements to be used on line 1
Transfinite Line{1} = numElementsX+1;
