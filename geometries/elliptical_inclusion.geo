// Radius of sphere
r_sphere = 100.0;

// Define points needed for a sphere
Point(1) = {0, 0, 0};
Point(2) = {r_sphere, 0, 0};
Point(3) = {0, r_sphere, 0};
Point(4) = {-r_sphere, 0, 0};
Point(5) = {0, -r_sphere, 0};
Point(6) = {0, 0, r_sphere};
Point(7) = {0, 0, -r_sphere};

// Circle in xy plane
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Circle in xz plane
Circle(5) = {6, 1, 4};
Circle(6) = {6, 1, 2};
Circle(7) = {7, 1, 4};
Circle(8) = {7, 1, 2};

// Circle in yz plane
Circle(9) = {6, 1, 3};
Circle(10) = {6, 1, 5};
Circle(11) = {7, 1, 3};
Circle(12) = {7, 1, 5};

// Create spherical surface
Line Loop(13) = {1, -9, 6};
Ruled Surface(14) = {13};

Line Loop(15) = {1, -11, 8};
Ruled Surface(16) = {15};

Line Loop(17) = {2, -5, 9};
Ruled Surface(18) = {17};

Line Loop(19) = {2, -7, 11};
Ruled Surface(20) = {19};

Line Loop(21) = {3, -10, 5};
Ruled Surface(22) = {21};

Line Loop(23) = {3, -12, 7};
Ruled Surface(24) = {23};

Line Loop(25) = {4, -6, 10};
Ruled Surface(26) = {25};

Line Loop(27) = {4, -8, 12};
Ruled Surface(28) = {27};

// Create spherical solid
Surface Loop(29) = {14, 16, 18, 20, 22, 24, 26, 28};
Volume(30) = {29};

