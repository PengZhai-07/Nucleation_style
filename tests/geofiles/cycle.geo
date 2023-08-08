// a geofile to generate structured transfinite quad mesh
//+
h1=80;
h2=800;
Point(1) = {0, 0, 0, h1};
Point(2) = {32000, 0, 0, h2};
Point(3) = {32000, 48000, 0, h2};
Point(4) = {0, 48000, 0, h1};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};

//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};

Physical Surface('mat1') = {1};

Physical Curve('bc1') = {1};
//+
Physical Curve('bc2') = {2};
//+
Physical Curve('bc3') = {3};
//+
Physical Curve('bc4') = {4};

// Transfinite Curve {1:4} = 3;
// Transfinite Curve {5,6,7} = 3;
// Transfinite Surface{1,2};

Recombine Surface{1};
Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2; // or 3
// RefineMesh;
