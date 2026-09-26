SetFactory("OpenCASCADE");

// outer radius and aspect ratio alpha = a/b, overridable with -setnumber
DefineConstant[ b = 1.0 ];
DefineConstant[ alpha = 0.5 ];

Circle(1) = {0, 0, 0, b};
Circle(2) = {0, 0, 0, alpha*b};
Curve Loop(1) = {1};
Curve Loop(2) = {2};
Plane Surface(1) = {1, 2};
