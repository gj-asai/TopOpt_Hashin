Point(1) = {0, 200, 0, 1.0};
Point(2) = {0, 100, 0, 1.0};
Point(3) = {0, 0, 0, 1.0};
Point(4) = {100, 0, 0, 1.0};
Point(5) = {200, 0, 0, 1.0};
Point(6) = {200, 100, 0, 1.0};
Point(7) = {100, 100, 0, 1.0};
Point(8) = {100, 200, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 7};
Line(10) = {4, 7};

Curve Loop(1) = {1, 9, 7, 8};
Curve Loop(2) = {2, 3, 10, -9};
Curve Loop(3) = {4, 5, 6, -10};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} = 25;
