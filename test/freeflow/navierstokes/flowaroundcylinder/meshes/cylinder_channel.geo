// channel measurements according to
// the DFG benchmark
//

// domain mesh size
cl_ = 0.01;
// cylinder surface mesh size
cl2_ = 0.002;

//+
Point(1) = {0.0, 0.0, 0, cl_};
//+
Point(2) = {0.0, 0.0, 0, cl_};
//+
Point(3) = {2.2, 0, 0, cl_};
//+
Point(4) = {0.0, 0.41, 0.0, cl_};
//+
Point(5) = {2.2, 0.41, 0.0, cl_};
//+
Point(6) = {0.15, 0.2, 0.0, cl2_};
//+
Point(7) = {0.25, 0.2, 0.0, cl2_};
//+
Point(8) = {0.2, 0.15, 0.0, cl2_};
//+
Point(9) = {0.2, 0.25, 0.0, cl2_};
//+
Point(10) = {0.2, 0.2, 0, cl2_};
//+
Line(1) = {4, 5};
//+
Line(2) = {5, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 4};
//+
Circle(5) = {6, 10, 9};
//+
Circle(6) = {9, 10, 7};
//+
Circle(7) = {7, 10, 8};
//+
Circle(8) = {8, 10, 6};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {6, 7, 8, 5};
//+
Plane Surface(1) = {1, 2};
