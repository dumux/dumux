cl_ = 1.0;
radius = 0.0002;
outer = 0.006;
box = 0.5*Sqrt(Pi*outer*outer);

//+
Point(1) = {0, 0, 0, cl_};
//+
Point(2) = {radius, 0, 0, cl_};
//+
Point(3) = {0, radius, 0, cl_};
//+
Circle(1) = {3, 1, 2};
//+
Point(4) = {box, 0, 0, cl_};
//+
Point(5) = {0, box, 0, cl_};
//+
Point(6) = {box, box, 0, cl_};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {5, 3};
//+
Curve Loop(1) = {5, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};

Field[1] = MathEval;
Field[1].F = Sprintf("(10-0.3)*%g/Log(%g/%g)*Log(Sqrt(x*x + y*y)/%g) + 0.3*%g", radius, outer, radius, radius, radius);
Background Field = 1;

