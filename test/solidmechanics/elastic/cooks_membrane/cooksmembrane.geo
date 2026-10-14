// Cook's membrane benchmark geometry
// Tapered panel with vertices at (0,0), (48,44), (48,60), (0,44)
// Left edge (x=0): clamped (Dirichlet zero)
// Right edge (x=48): shear load (Neumann)
// Top/bottom edges: traction-free (zero Neumann)

h = 4.0;

Point(1) = {0,  0,  0, h};
Point(2) = {48, 44, 0, h};
Point(3) = {48, 60, 0, h};
Point(4) = {0,  44, 0, h};

Line(1) = {1, 2}; // bottom edge
Line(2) = {2, 3}; // right edge (shear load)
Line(3) = {3, 4}; // top edge
Line(4) = {4, 1}; // left edge (clamped)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
