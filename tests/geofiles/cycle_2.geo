// a geofile to generate structured transfinite quad mesh
//+
h1=80;
h2=800;
Point(1) = {0, 0, 0, h1};
Point(2) = {8000, 0, 0, h1};
Point(3) = {8000, 48000, 0, h1};
Point(4) = {0, 48000, 0, h1};
Point(5) = {32000, 0, 0, h2};
Point(6) = {32000, 48000, 0, h2};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

//+
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {-5, 2, -7,-6};
//+
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Surface('mat1') = {1};
Physical Surface('mat2') = {2};

Physical Curve('bc1') = {1,5};
//+
Physical Curve('bc2') = {6};
//+
Physical Curve('bc3') = {3,7};
//+
Physical Curve('bc4') = {4};

// Transfinite Curve {1:4} = 3;
// Transfinite Curve {5,6,7} = 3;
// Transfinite Surface{1,2};

Recombine Surface{1,2};
Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2; // or 3
// RefineMesh;
//+

