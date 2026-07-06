// Hysing et al. (2009) FeatFlow rising-bubble benchmark, HALF domain [0,0.5] x [0,2].
// Case 1 is left-right symmetric, so we simulate only the left half and place a
// symmetry plane at x=0.5 (the bubble centre). The benchmark's free-slip side-wall BC
// (no penetration u_x=0 + zero tangential shear) IS the symmetry condition, so the
// x=0.5 boundary needs no special treatment - it is handled as a normal "side wall".
// This enforces exact symmetry by construction (no mesh-seeded drift) and halves the cost.
// The no-slip floor is at y=0; the bubble is centred at (0.5,0.5), r=0.25, so its left
// half occupies x in [0.25,0.5] against the symmetry plane.
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthMax = 0.1; // base h ~ 0.1 (AMR refines the interface)
Rectangle(2) = {0, 0, 0, 0.5, 2, 0};
