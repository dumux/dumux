a = 1;             // width of the domain
b = 1;             // height of the domain

numCells = 50;
dx_m = (a+b)/2.0/numCells; // discretization length in the matrix
dx_f = dx_m*0.5;           // discretization length in the fracture


// Points
Point(1) = {0.0,   0.0, 0.0, dx_m};
Point(2) = {a,     0.0, 0.0, dx_m};
Point(3) = {a,     b,   0.0, dx_m};
Point(4) = {0.0,   b,   0.0, dx_m};

Point(5) = {0.25, 0.5*b, 0.0, dx_f};
Point(6) = {0.75, 0.5*b, 0.0, dx_f};

// domain outline
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(1) = {1:4};

// domain surface
Plane Surface(1) = {1};

// fracture line
Line(6) = {5, 6};
Line{6} In Surface{1};

// make fracture & surface physical entities
Physical Surface(1) = {1};
Physical Line(1) = {6};
