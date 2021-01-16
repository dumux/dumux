lc = 1.0;
lc2 = 0.25;
Point(1) = {90.0, 46.409, 0.0, lc2};
Point(2) = {93.167, 45.924, 0.0, lc2};
Point(3) = {96.876, 45.356, 0.0, lc2};
Point(4) = {98.783, 45.064, 0.0, lc2};
Point(5) = {95.088, 40.0, 0.0, lc2};
Point(6) = {95.962, 42.617, 0.0, lc2};
Point(7) = {97.48, 47.168, 0.0, lc2};
Point(8) = {93.626, 50.0, 0.0, lc2};
Point(9) = {93.043, 44.818, 0.0, lc};
Point(10) = {99.3, 43.493, 0.0, lc};
Point(11) = {93.422, 41.941, 0.0, lc};
Point(12) = {90.0, 40.0, 0.0, lc};
Point(13) = {90.0, 50.0, 0.0, lc};
Point(14) = {100.0, 40.0, 0.0, lc};
Point(15) = {100.0, 50.0, 0.0, lc};

// for some reason the coordinates are not beginning in the origin
// we translate all points to the origin (Point 12 is the lower left corner)
Translate {-90, -40, 0} { Point{1:15}; }

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {5, 6};
Line(5) = {6, 3};
Line(6) = {3, 7};
Line(7) = {8, 2};
Line(8) = {2, 9};
Line(9) = {10, 6};
Line(10) = {6, 11};

Line(11) = {12, 1};
Line(12) = {1, 13};
Line(13) = {13, 8};
Line(14) = {8, 15};
Line(15) = {15, 14};
Line(16) = {14, 5};
Line(17) = {5, 12};
Line Loop(1) = {11, 12, 13, 14, 15, 16, 17};
Plane Surface(1) = {1};

Line{1} In Surface{1};
Line{2} In Surface{1};
Line{3} In Surface{1};
Line{4} In Surface{1};
Line{5} In Surface{1};
Line{6} In Surface{1};
Line{7} In Surface{1};
Line{8} In Surface{1};
Line{9} In Surface{1};
Line{10} In Surface{1};
Line{11} In Surface{1};
Line{12} In Surface{1};
Line{13} In Surface{1};
Line{14} In Surface{1};

Physical Line(1) = {1:10};
Physical Surface(1) = {1};
