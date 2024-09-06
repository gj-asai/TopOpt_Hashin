// Gmsh project created on Sun Aug 11 23:06:30 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {60, 0, 0, 1.0};
//+
Point(3) = {60, 20, 0, 1.0};
//+
Point(4) = {0, 20, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 3} = 61 Using Progression 1;
//+
Transfinite Curve {2, 4} = 21 Using Progression 1;
