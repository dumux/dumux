a = 500;        // width of the domain
b = 500;        // height of the domain
alpha = 20.0;   // fracture inclination angle [degrees]
df = 10;        // length of the fracture

numCells = 30;
numCellsOnFracture = 30;

dx_m = (a+b)/2.0/numCells;            // discretization length in the matrix
dx_f = df/numCellsOnFracture;         // discretization length in the fracture

// Domain boundary points
Point(1) = {-a/2.0, -b/2.0, 0.0, dx_m};
Point(2) = { a/2.0, -b/2.0, 0.0, dx_m};
Point(3) = { a/2.0,  b/2.0, 0.0, dx_m};
Point(4) = {-a/2.0,  b/2.0, 0.0, dx_m};

// Fracture corners
dx = df*Cos(alpha*Pi/180.0);
dy = df*Sin(alpha*Pi/180.0);
Point(5) = {-dx*0.5,  dy*0.5, 0.0, dx_f};
Point(6) = { dx*0.5, -dy*0.5, 0.0, dx_f};

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
