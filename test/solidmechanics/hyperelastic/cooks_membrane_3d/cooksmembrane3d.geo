// 3D Cook's membrane benchmark geometry
// Tapered panel extruded in z: cross-section (0,0)-(48,44)-(48,60)-(0,44), depth=1 mm
//
// BCs:
//   Left face  (x=0):  clamped (u=0)
//   Right face (x=48): shear load T=(0,P₀,0)
//   z=0 face:          symmetry (u_z=0)  → plane-strain equivalent
//   All other faces:   traction-free
//
// Same cross-section as the 2D test; depth = 1 (same units = mm).
// The plane-strain z=0 BC makes this directly comparable to the 2D result
// and to the SPP-1748 3D p-FEM reference (Table 2.9).

h = 4.0;
depth = 1.0;

// 2D cross-section vertices at z=0
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

// Extrude the cross-section by `depth` in z
out[] = Extrude {0, 0, depth} {
  Surface{1};
};
// out[0] = top surface (z=depth)
// out[1] = volume
// out[2..5] = extruded side surfaces (from lines 1..4 in order)

// Physical groups for boundary identification
Physical Volume("volume") = {out[1]};
Physical Surface("bottom_zface") = {1};          // z=0
Physical Surface("top_zface")    = {out[0]};     // z=depth
Physical Surface("bottom_edge")  = {out[2]};     // from Line 1 (y≈0, tilted)
Physical Surface("right_face")   = {out[3]};     // from Line 2 (x=48, shear load)
Physical Surface("top_edge")     = {out[4]};     // from Line 3 (y≈60, traction-free)
Physical Surface("left_face")    = {out[5]};     // from Line 4 (x=0, clamped)
