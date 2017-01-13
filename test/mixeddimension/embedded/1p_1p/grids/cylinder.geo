size = 0.1;
Point(1) = {0, 0, 0, size};
Point(2) = {1, 0, 0, size};
Point(3) = {0, 1, 0, size};
Point(4) = {0, -1, 0, size};
Point(5) = {-1, 0, 0, size};
Point(6) = {-1, 0, 1, size};
Point(7) = {0, 0, 1, size};
Point(8) = {0, 1, 1, size};
Point(13) = {1, 0, 1, size};
Point(18) = {0, -1, 1, size};

Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 3};
Circle(8) = {6, 7, 8};
Circle(9) = {8, 7, 13};
Circle(10) = {13, 7, 18};
Circle(11) = {18, 7, 6};

Line(13) = {5, 6};
Line(14) = {3, 8};
Line(18) = {2, 13};
Line(22) = {4, 18};

Line Loop(6) = {4, 1, 2, 3};
Ruled Surface(6) = {6};
Line Loop(15) = {4, 14, -8, -13};
Ruled Surface(15) = {15};
Line Loop(19) = {1, 18, -9, -14};
Ruled Surface(19) = {19};
Line Loop(23) = {2, 22, -10, -18};
Ruled Surface(23) = {23};
Line Loop(27) = {3, 13, -11, -22};
Ruled Surface(27) = {27};
Line Loop(28) = {8, 9, 10, 11};
Ruled Surface(28) = {28};
Surface Loop(1) = {6, 28, 15, 19, 23, 27};

Volume(1) = {1};
