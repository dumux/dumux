///////////////////////////////////////////////////////
// user input /////////////////////////////////////////
///////////////////////////////////////////////////////

//fault
x_fault_bottom = 400;
y_fault_bottom = 500;
x_fault_top = 600;
y_fault_top = 1500;
faultwidth = 12.5; // half of the distance is added to the left and right
y_tolerance = 5; // to use only quadrilateral elements (gets rid of triangle)

//domain
x_min = 0.0;
x_max = 2000;
y_min = 0.0;
y_max = 2000;

//cells in x and y direction
x_cells_faultwidth = 5;
delta_x = 100;
delta_y_underburden = 100;
delta_y_overburden = 100;
delta_y_fault = 0.5;

// cells_z = 1;
// depth = 1; //in z-direction

///////////////////////////////////////////////////////
anglefault = Atan( (y_fault_top - y_fault_bottom) / (x_fault_top - x_fault_bottom) );

correction = faultwidth * Cos(anglefault);

x_correction = correction * Cos(anglefault);
y_correction = correction * Sin(anglefault);

x_tolerance = y_tolerance / Tan(anglefault);

Point(1)  = {x_min,                                                             y_fault_bottom,                                 0, 1.0};
Point(11) = {x_min,                                                             y_fault_bottom + y_correction + y_tolerance,    0, 1.0};
Point(2)  = {x_fault_bottom-faultwidth/2,                                       y_fault_bottom,                                 0, 1.0};
Point(21) = {x_fault_bottom -faultwidth/2    + x_correction + x_tolerance,      y_fault_bottom + y_correction + y_tolerance,    0, 1.0};
Point(3)  = {x_fault_bottom+faultwidth/2,                                       y_fault_bottom,                                 0, 1.0};
Point(31) = {x_fault_bottom +faultwidth/2    + x_tolerance,                     y_fault_bottom                + y_tolerance,    0, 1.0};
Point(4)  = {x_max,                                                             y_fault_bottom,                                 0, 1.0};
Point(41) = {x_max,                                                             y_fault_bottom                + y_tolerance,    0, 1.0};

Point(5)  = {x_min,                                                             y_fault_top,                                    0, 1.0};
Point(51) = {x_min,                                                             y_fault_top                   - y_tolerance,    0, 1.0};
Point(6)  = {x_fault_top-faultwidth/2,                                          y_fault_top,                                    0, 1.0};
Point(61) = {x_fault_top    -faultwidth/2    - x_tolerance,                     y_fault_top                   - y_tolerance,    0, 1.0};
Point(7)  = {x_fault_top+faultwidth/2,                                          y_fault_top,                                    0, 1.0};
Point(71) = {x_fault_top    +faultwidth/2    - x_correction - x_tolerance,      y_fault_top    - y_correction - y_tolerance,    0, 1.0};
Point(8)  = {x_max,                                                             y_fault_top,                                    0, 1.0};
Point(81) = {x_max,                                                             y_fault_top    - y_correction - y_tolerance,    0, 1.0};

Line(101) = {1, 2};
Line(102) = {2, 3};
Line(103) = {3, 4};

Line(104) = {5, 6};
Line(105) = {6, 7};
Line(106) = {7, 8};

Line(201) = {11, 51};
Line(202) = {21, 61};
Line(203) = {31, 71};
Line(204) = {41, 81};

Line(1011) = {11, 21};
Line(1021) = {21, 31};
Line(1031) = {31, 41};
Line(1041) = {51, 61};
Line(1051) = {61, 71};
Line(1061) = {71, 81};

Line(2011) = {1, 11};
Line(2012) = {51, 5};
Line(2021) = {2, 21};
Line(2022) = {61, 6};
Line(2031) = {3, 31};
Line(2032) = {71, 7};
Line(2041) = {4, 41};
Line(2042) = {81, 8};

Line Loop(1) = {1011, 202, -1041, -201};
Plane Surface(1) = {1};
Line Loop(11) = {101, 2021, -1011, -2011};
Plane Surface(11) = {11};
Line Loop(12) = {1041, 2022, -104, -2012};
Plane Surface(12) = {12};

Line Loop(2) = {1021, 203, -1051, -202};
Plane Surface(2) = {2};
Line Loop(21) = {102, 2031, -1021, -2021};
Plane Surface(21) = {21};
Line Loop(22) = {1051, 2032, -105, -2022};
Plane Surface(22) = {22};

Line Loop(3) = {1031, 204, -1061, -203};
Plane Surface(3) = {3};
Line Loop(31) = {103, 2041, -1031, -2031};
Plane Surface(31) = {31};
Line Loop(32) = {1061, 2042, -106, -2032};
Plane Surface(32) = {32};

//bottom
x_cells_leftoffault  = (x_fault_bottom - x_min) / delta_x + 1;
x_cells_rightoffault = (x_max - (x_fault_bottom + faultwidth)) / delta_x + 1;

y_cells_fault = (y_fault_top - y_fault_bottom) / delta_y_fault + 1;

Transfinite Line {101} = x_cells_leftoffault;
Transfinite Line {1011} = x_cells_leftoffault;
Transfinite Line {102} = x_cells_faultwidth + 1;
Transfinite Line {1021} = x_cells_faultwidth + 1;
Transfinite Line {103} = x_cells_rightoffault;
Transfinite Line {1031} = x_cells_rightoffault;

//top
Transfinite Line {104} = x_cells_leftoffault;
Transfinite Line {1041} = x_cells_leftoffault;
Transfinite Line {105} = x_cells_faultwidth + 1;
Transfinite Line {1051} = x_cells_faultwidth + 1;
Transfinite Line {106} = x_cells_rightoffault;
Transfinite Line {1061} = x_cells_rightoffault;

//fault
// fault_cells = (fault_top_y-fault_bottom_y)/delta_y*cells_y;
Transfinite Line {201} = y_cells_fault;
Transfinite Line {2011} = 2;
Transfinite Line {2012} = 2;
Transfinite Line {202} = y_cells_fault;
Transfinite Line {2021} = 2;
Transfinite Line {2022} = 2;
Transfinite Line {203} = y_cells_fault;
Transfinite Line {2031} = 2;
Transfinite Line {2032} = 2;
Transfinite Line {204} = y_cells_fault;
Transfinite Line {2041} = 2;
Transfinite Line {2042} = 2;

Transfinite Surface {1, 11, 12, 2, 21, 22, 3, 31, 32};
Recombine Surface {1, 11, 12, 2, 21, 22, 3, 31, 32};

cells_y_overburden  = (y_max - y_fault_top)     / delta_y_overburden + 1;
cells_y_underburden = (y_fault_bottom - y_min)  / delta_y_underburden +1;

//Extrude top and bottom domain
s[] = Extrude{0,y_min-y_fault_bottom,0}{Line{101,102,103}; Layers{cells_y_underburden}; Recombine;};
s2[] = Extrude{0,y_max-y_fault_top,0}{Line{104,105,106}; Layers{cells_y_overburden}; Recombine;};

//Extrude one cell in depth
// Extrude{0,0,depth}{Surface{s[1],s[5],s2[1],s2[5],2,3}; Layers{cells_z}; Recombine;};

