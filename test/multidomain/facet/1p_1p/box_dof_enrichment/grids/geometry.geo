// domain measures
xmax = 1.0;
ymax = 3.0;
inlet_size = 0.2;
mesh_size = 0.05;

Point(0) = {0, 0, 0, mesh_size};
Point(1) = {xmax, 0, 0, mesh_size};

Point(2) = {0, ymax/2, 0, mesh_size};
Point(3) = {inlet_size, ymax/2, 0, mesh_size};
Point(4) = {xmax - inlet_size, ymax/2, 0, mesh_size};
Point(5) = {xmax, ymax/2, 0, mesh_size};

Point(6) = {inlet_size, ymax, 0, mesh_size};
Point(7) = {xmax - inlet_size, ymax, 0, mesh_size};

Line(0) = {0, 1};
Line(1) = {0, 2};

Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};

Line(5) = {3, 6};
Line(6) = {4, 7};

Line(7) = {6, 7};
Line(8) = {1, 5};

Curve Loop(1) = {0, 8, -4, -3, -2, -1};
Plane Surface(1) = {1};

Curve Loop(2) = {3, 6, -7, -5};
Plane Surface(2) = {2};

// to not extend the fracture beyond the top-domain boundary,
// you may use the command-line argument "-setnumber extend_fracture 0"
DefineConstant[ extend_fracture = 1 ];

Physical Line(1) = {3};
If (extend_fracture)
    Physical Line(2) = {2, 4};
EndIf

Physical Surface(1) = {1, 2};
