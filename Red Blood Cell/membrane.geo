s=2.0;
r=10;

Point(1) = {0, 0, 0, s};
Point(2) = {r, 0, 0, s};
Point(3) = {-r, 0, 0, s};
Point(4) = {0, r, 0, s};
Point(5) = {0, -r, 0, s};
Point(6) = {0, 0, r, s};
Point(7) = {0, 0, -r, s};//+
Circle(1) = {3, 1, 6};
//+
Circle(2) = {6, 1, 2};
//+
Circle(3) = {2, 1, 7};
//+
Circle(4) = {7, 1, 3};
//+
Circle(5) = {4, 1, 7};
//+
Circle(6) = {7, 1, 5};
//+
Circle(7) = {5, 1, 6};
//+
Circle(8) = {6, 1, 4};
//+
Circle(9) = {3, 1, 4};
//+
Circle(10) = {4, 1, 2};
//+
Circle(11) = {2, 1, 5};
//+
Circle(12) = {5, 1, 3};//+
Line Loop(13) = {9, -8, -1};
//+
Ruled Surface(14) = {13};
//+
Line Loop(15) = {8, 10, -2};
//+
Ruled Surface(16) = {15};
//+
Line Loop(17) = {10, 3, -5};
//+
Ruled Surface(18) = {17};
//+
Line Loop(19) = {5, 4, 9};
//+
Ruled Surface(20) = {19};
//+
Line Loop(21) = {6, 12, -4};
//+
Ruled Surface(22) = {21};
//+
Line Loop(23) = {12, 1, -7};
//+
Ruled Surface(24) = {23};
//+
Line Loop(25) = {2, 11, 7};
//+
Ruled Surface(26) = {25};
//+
Line Loop(27) = {11, -6, -3};
//+
Ruled Surface(28) = {27};
//+
Surface Loop(29) = {22, 28, 26, 16, 14, 20, 18, 24};
//+
Volume(30) = {29};
//+
Physical Surface(31) = {28, 26, 16, 18, 20, 22, 24, 14};
//+
Physical Volume(32) = {30};
