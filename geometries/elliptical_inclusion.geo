// Radius of sphere
r_sphere = 100.0;

// Major axis of ellipsoid
a = 10.0;
b = 5.0;
c = 2.0;

// Define mesh control parameters
lc_inc = 0.1; // Lattice control for inclusion
lc_sph = 5.0; // Lattice control for outer sphere

// Define points needed for a sphere
Point(1) = {0, 0, 0, lc_inc};
Point(2) = {r_sphere, 0, 0, lc_sph};
Point(3) = {0, r_sphere, 0, lc_sph};
Point(4) = {-r_sphere, 0, 0, lc_sph};
Point(5) = {0, -r_sphere, 0, lc_sph};
Point(6) = {0, 0, r_sphere, lc_sph};
Point(7) = {0, 0, -r_sphere, lc_sph};

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

// Define points for ellipsoid start number on 100
Point(101) = {a, 0, 0, lc_inc};
Point(102) = {-a, 0, 0, lc_inc};
Point(103) = {0, b, 0, lc_inc};
Point(104) = {0, -b, 0, lc_inc};
Point(105) = {0, 0, c, lc_inc};
Point(106) = {0, 0, -c, lc_inc};

// Create ellipse in xy plane
Ellipsis(107) = {101, 1, 102, 103};
Ellipsis(108) = {102, 1, 101, 104};
Ellipsis(109) = {101, 1, 102, 104};
Ellipsis(110) = {102, 1, 101, 103};

// Create ellipse in the xz plane
Ellipsis(111) = {101, 1, 102, 105};
Ellipsis(112) = {101, 1, 102, 106};
Ellipsis(113) = {102, 1, 101, 105};
Ellipsis(114) = {102, 1, 101, 106};

// Create ellipse in the yz plane
Ellipsis(115) = {103, 1, 104, 105};
Ellipsis(116) = {103, 1, 104, 106};
Ellipsis(117) = {104, 1, 103, 105};
Ellipsis(118) = {104, 1, 103, 106};

// Create surface of the ellipsoid
Line Loop(119) = {107, 115, -111};
Ruled Surface(120) = {119};

Line Loop(121) = {107, 116, -112};
Ruled Surface(122) = {121};

Line Loop(123) = {108, 117, -113};
Ruled Surface(124) = {123};

Line Loop(125) = {108, 118, -114};
Ruled Surface(126) = {125};

Line Loop(127) = {109, 117, -111};
Ruled Surface(128) = {127};

Line Loop(129) = {109, 118, -112};
Ruled Surface(130) = {129};

Line Loop(131) = {110, 115, -113};
Ruled Surface(132) = {131};

Line Loop(133) = {110, 116, -114};
Ruled Surface(134) = {133};
 
// Create spherical solid
Surface Loop(29) = {14, 16, 18, 20, 22, 24, 26, 28};

// Create elliptical surface
Surface Loop(135) = {120, 122, 124, 126, 128, 130, 132, 134};

Volume(30) = {29, 135};
Volume(136) = {135};

// Add physical entities
Physical Volume("matrix") = {30};
Physical Volume("inclusion") = {136};
Physical Surface("outer") = {29};

