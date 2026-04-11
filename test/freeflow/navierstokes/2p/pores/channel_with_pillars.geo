// ============================================================
//  2D GMSH geometry: Gas channel + pillar box + cone exit
// ============================================================
//
//  Layout (not to scale):
//
//   |<------------ channel_L ------------->|
//   +--------------------------------------+  ^
//   |          GAS CHANNEL                 |  | channel_H
//   |    |<----- box_W ------>|            |  v
//   +----+--------------------+------------+
//        |   PILLAR  BOX      |
//        |  o  o  o  o  o     |
//        |  o  o  o  o  o     |  box_H
//        |  o  o  o  o  o     |
//        +---\           /----+
//              \         /      <- cone walls
//               +-[exit]-+       exit width = 0.1 * box_W
//
// ============================================================

// Use OpenCASCADE kernel for boolean operations and primitives
SetFactory("OpenCASCADE");

// ---------- USER PARAMETERS (edit these) ----------

channel_L   = 1.0;    // Total length of the gas channel
channel_H   = 0.1;    // Height of the gas channel

box_W       = 0.5 * channel_L;  // Width of the pillar box  (= 50% of channel_L)
box_H       = 0.4;    // Height of the pillar box

pillar_R    = 0.025;  // Radius of each circular pillar
pillar_nx   = 6;      // Number of pillars in x-direction
pillar_ny   = 4;      // Number of pillars in y-direction

// Cone / nozzle at the bottom of the box
exit_W      = 0.1 * box_W;   // Width of the flat exit slot
y_offset_exit = 0.1 * box_H;

// Mesh sizes
lc_channel  = 0.03;  // Mesh size in the open channel regions
lc_box      = 0.02;  // Mesh size inside the pillar box
lc_pillar   = 0.001; // Mesh size around pillars
lc_exit     = 0.005; // Mesh size at the cone exit

// Distance from first row of pillars to gas channel
pillar_margin = 0.1;

// --------------------------------------------------
// Derived coordinates
// --------------------------------------------------

// Channel baseline (y = 0 is the bottom of the channel)
x_chan_left  = 0.0;
x_chan_right = channel_L;
y_chan_bot   = 0.0;
y_chan_top   = channel_H;

// Box is centred in x on the channel
x_box_left   = 0.5 * (channel_L - box_W);
x_box_right  = x_box_left + box_W;
y_box_top    = y_chan_bot;          // Box sits below the channel
y_box_bot    = y_chan_bot - box_H;

// Cone exit (centred on the box)
x_exit_left  = 0.5 * (x_box_left + x_box_right) - 0.5 * exit_W;
x_exit_right = x_exit_left + exit_W;
y_exit       = y_box_bot - y_offset_exit;   // bottom of the box / top of cone
// Cone apex level is the same as box bottom; the exit IS the flat bottom

// --------------------------------------------------
// Points – channel outer rectangle
// --------------------------------------------------
Point(1)  = {x_chan_left,  y_chan_bot, 0, lc_channel};
Point(2)  = {x_box_left,   y_chan_bot, 0, lc_box};
Point(3)  = {x_box_right,  y_chan_bot, 0, lc_box};
Point(4)  = {x_chan_right,  y_chan_bot, 0, lc_channel};
Point(5)  = {x_chan_right,  y_chan_top, 0, lc_channel};
Point(6)  = {x_chan_left,   y_chan_top, 0, lc_channel};

// Points – pillar box sides (vertical walls hang down from channel bottom)
Point(7)  = {x_box_left,   y_box_bot,  0, lc_box};
Point(8)  = {x_box_right,  y_box_bot,  0, lc_box};

// Points – cone (the box bottom has a trapezoidal cut)
//   Left cone foot  -> exit left
//   Right cone foot -> exit right
Point(9)  = {x_exit_left,  y_exit,     0, lc_exit};
Point(10) = {x_exit_right, y_exit,     0, lc_exit};

// --------------------------------------------------
// Lines – channel top + sides
// --------------------------------------------------
Line(1)  = {1, 2};   // channel bottom, left of box
Line(2)  = {3, 4};   // channel bottom, right of box
Line(3)  = {4, 5};   // channel right wall
Line(4)  = {5, 6};   // channel top wall
Line(5)  = {6, 1};   // channel left wall

// Lines – interface between channel and box (shared boundary)
Line(6)  = {2, 7};   // box left wall
Line(7)  = {7, 9};   // box bottom, left of exit  (cone left wall)
Line(8)  = {9, 10};  // exit slot (outflow boundary)
Line(9)  = {10, 8};  // box bottom, right of exit (cone right wall)
Line(10) = {8, 3};   // box right wall

// --------------------------------------------------
// Curve Loops & Surfaces
// --------------------------------------------------

// --- Channel (without the box opening) ---
// Bottom of channel has a gap where the box connects (points 2..3)
// We close it with a virtual line across the box top
Line(11) = {2, 3};   // virtual top of box / bottom of channel opening

Curve Loop(1) = {1, 11, 2, 3, 4, 5};
Plane Surface(1) = {1};          // GAS CHANNEL surface

// --- Pillar box interior ---
// Boundary: left wall (6), cone-left (7), exit (8), cone-right (9),
//           right wall (10), and back up through the channel opening (-11)
Curve Loop(2) = {6, 7, 8, 9, 10, -11};
Plane Surface(2) = {2};          // BOX surface (pillars will be subtracted)

// --------------------------------------------------
// Circular pillars (regular grid inside the box)
// --------------------------------------------------

// Usable interior of the box (keep pillars away from walls)
x_inner_left  = x_box_left  + pillar_R;
x_inner_right = x_box_right - pillar_R;
y_inner_top   = y_box_top   - pillar_R - pillar_margin;  // with margin from gas channel
y_inner_bot   = y_box_bot   + pillar_R;  // above the cone region

// Spacing: diameter + gap of one radius (3 * pillar_R)
dx_step = 3 * pillar_R;
// Hexagonal packing for space-filling alternating rows
dy_step = Sqrt(3) / 2 * dx_step;

// Space needed for shifted rows (extends by 0.5*dx_step)
width_needed = (pillar_nx - 1) * dx_step + 0.5 * dx_step;
x_offset = (x_inner_right - x_inner_left - width_needed) / 2;

height_needed = (pillar_ny - 1) * dy_step;
y_offset = (y_inner_top - y_inner_bot - height_needed) / 2;

pillar_disks[] = {};    // will collect all pillar disk tags
pillar_centers_x[] = {};  // store pillar center x coordinates
pillar_centers_y[] = {};  // store pillar center y coordinates

base_sf  = 10;    // starting surface index for pillars

For j In {0 : pillar_ny - 1}
  For i In {0 : pillar_nx - 1}

    // Alternating shift: odd rows shifted by half spacing
    cx = x_inner_left + x_offset + i * dx_step + (j % 2) * 0.5 * dx_step;
    cy = y_inner_bot  + y_offset + j * dy_step;

    idx = j * pillar_nx + i;   // 0-based pillar index

    // Store pillar center coordinates
    pillar_centers_x[] += cx;
    pillar_centers_y[] += cy;

    // Create disk using OpenCASCADE primitive
    Disk(base_sf + idx) = {cx, cy, 0, pillar_R, pillar_R};

    // Collect disk tags for box subtraction
    pillar_disks[] += (base_sf + idx);

  EndFor
EndFor

// Subtract all pillar disks from the box surface using BooleanDifference
If (#pillar_disks[] > 0)
  BooleanDifference{ Surface{2}; Delete; }{ Surface{pillar_disks[]}; Delete; }
EndIf

// Create mesh refinement fields for each pillar
// Each field has linear ramp: fine at center, coarse at 2*R distance
ramp_distance = 2 * pillar_R;

// Create temporary lines through pillar centers for Attractor fields
field_list[] = {};
For idx In {0 : #pillar_centers_x[] - 1}
  pt_id = 6000 + idx;
  line_id = 7000 + idx;

  // Create a point at pillar center
  Point(pt_id) = {pillar_centers_x[idx], pillar_centers_y[idx], 0, 1.0};

  // Create a small line segment as anchor for Attractor field
  Point(pt_id + 100) = {pillar_centers_x[idx] + 0.001, pillar_centers_y[idx], 0, 1.0};
  Line(line_id) = {pt_id, pt_id + 100};

  field_id = idx + 1;

  // Attractor field to pillar center
  Field[field_id] = Attractor;
  Field[field_id].EdgesList = {line_id};
  Field[field_id].NumPointsPerCurve = 20;

  // Threshold field for linear ramp
  thresh_id = 1000 + idx;
  Field[thresh_id] = Threshold;
  Field[thresh_id].InField = field_id;
  Field[thresh_id].SizeMin = lc_pillar;
  Field[thresh_id].SizeMax = lc_box;
  Field[thresh_id].DistMin = 0;
  Field[thresh_id].DistMax = ramp_distance;

  field_list[] += {thresh_id};
EndFor

// Combine all pillar fields with minimum operation
If (#field_list[] > 0)
  Field[2000] = Min;
  For i In {0 : #field_list[] - 1}
    Field[2000].FieldsList += field_list[i];
  EndFor
  Background Field = 2000;
EndIf

// --------------------------------------------------
// Physical groups (for boundary conditions / export)
// --------------------------------------------------

// Boundaries of the channel
Physical Curve("inlet",         10001) = {5};      // left wall of channel
Physical Curve("outlet_chan",   10002) = {3};      // right wall of channel
Physical Curve("top_wall",      10003) = {4};      // top wall
Physical Curve("chan_bot_left", 10004) = {1};      // channel bottom, left
Physical Curve("chan_bot_right",10005) = {2};      // channel bottom, right

// Surfaces
Physical Surface("gas_channel", 20001) = {1};
Physical Surface("pillar_box",  20002) = {2};

// --------------------------------------------------
// Meshing options
// --------------------------------------------------
Mesh.Algorithm       = 6;   // Frontal-Delaunay (good for circles)
Mesh.RecombineAll    = 0;   // triangles  (set to 1 for quads)
Mesh.CharacteristicLengthMin = lc_pillar;
Mesh.CharacteristicLengthMax = lc_channel;

// Generate 2-D mesh
// Mesh 2;
