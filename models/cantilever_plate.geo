Point(1) = {0, 0, 0, 1.0};
Point(2) = {70, 0, 0, 1.0};
Point(3) = {70, 70, 0, 1.0};
Point(4) = {0, 70, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {4, 1, 2, 3};

Plane Surface(1) = {1};
Transfinite Surface {1};
Transfinite Curve {2, 1, 4, 3} = 31 Using Progression 1;
