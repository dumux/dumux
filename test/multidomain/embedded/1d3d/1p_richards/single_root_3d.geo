radius = 0.0003;
outer = 0.006;
box = 0.5*Sqrt(Pi*(outer*outer));
height = 0.001;

//+
SetFactory("OpenCASCADE");

//+
Box(1) = {-box, -box, -height, 2*box, 2*box, height};

Field[1] = MathEval;
Field[1].F = Sprintf("(2-0.1)*%g/Log(%g/%g)*Log(Sqrt(x*x + y*y)/%g) + 0.1*%g", radius, outer, radius, radius, radius);

Field[2] = MathEval;
Field[2].F = Sprintf("0.1*%g", radius);

Field[3] = Max;
Field[3].FieldsList = {1, 2};

Background Field = 3;
