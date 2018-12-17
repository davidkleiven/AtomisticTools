// Radius of sphere
L_box = 30.0;

// Major axis of ellipsoid
a = 10.0;
b = 5.0;
c = 2.0;

// Define mesh control parameters
lc_inc = 0.1; // Lattice control for inclusion
lc_box = 5.0; // Lattice control for outer sphere

// Define point at origin
Point(1) = {0, 0, 0, lc_inc};


Point(2) = {-L_box/2, -L_box/2, -L_box/2, lc_box};
Point(3) = {L_box/2, -L_box/2, -L_box/2, lc_box};
Point(4) = {-L_box/2, L_box/2, -L_box/2, lc_box};
Point(5) = {L_box/2, L_box/2, -L_box/2, lc_box};
Point(6) = {-L_box/2, -L_box/2, L_box/2, lc_box};
Point(7) = {L_box/2, -L_box/2, L_box/2, lc_box};
Point(8) = {-L_box/2, L_box/2, L_box/2, lc_box};
Point(9) = {L_box/2, L_box/2, L_box/2, lc_box};

// Define connecting lines
Line(1) = {2, 3};
Line(2) = {2, 6};
Line(3) = {2, 4};
Line(4) = {3, 5};
Line(5) = {3, 7};
Line(6) = {4, 5};
Line(7) = {4, 8};
Line(8) = {5, 9};
Line(9) = {6, 8};
Line(10) = {6, 7};
Line(11) = {7, 9};
Line(12) = {8, 9};

// Define planes
Line Loop(134) = {9, -7, -3, 2};
Plane Surface(135) = {134};
Line Loop(135) = {12, -11, -10, 9};
Plane Surface(136) = {135};
Line Loop(136) = {7, 12, -8, -6};
Plane Surface(137) = {136};
Line Loop(137) = {3, 6, -4, -1};
Plane Surface(138) = {137};
Line Loop(138) = {4, 8, -11, -5};
Plane Surface(139) = {138};
Line Loop(139) = {10, -5, -1, 2};
Plane Surface(140) = {139};

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
Surface(120) = {119};

Line Loop(121) = {107, 116, -112};
Surface(122) = {121};

Line Loop(123) = {108, 117, -113};
Surface(124) = {123};

Line Loop(125) = {108, 118, -114};
Surface(126) = {125};

Line Loop(127) = {109, 117, -111};
Surface(128) = {127};

Line Loop(129) = {109, 118, -112};
Surface(130) = {129};

Line Loop(131) = {110, 115, -113};
Surface(132) = {131};

Line Loop(133) = {110, 116, -114};
Surface(134) = {133};

// Create elliptical surface
Surface Loop(135) = {120, 122, 124, 126, 128, 130, 132, 134};

Surface Loop(200) = {135, 136, 137, 139, 138, 140};
Volume(137) = {135, 200};
Volume(136) = {135};

// Create periodic surfaces
Periodic Surface{135} = {139} Translate(-L_box, 0, 0);
Periodic Surface {137} = {140} Translate(0, L_box, 0); // y-direction
Periodic Surface {136} = {138} Translate(0, 0, L_box); // z-direction

// Add physical entities
Physical Volume("matrix") = {137};
Physical Volume("inclusion") = {136};
Physical Surface("outer") = {135, 136, 137, 139, 138, 140};

