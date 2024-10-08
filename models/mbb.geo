Point(1) = {0, 0, 0, 1.0};
Point(2) = {168, 0, 0, 1.0};
Point(3) = {168, 80, 0, 1.0};
Point(4) = {0, 80, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};

Transfinite Surface {1};
Transfinite Curve {1, 3} = 85 Using Progression 1;
Transfinite Curve {2, 4} = 41 Using Progression 1;
