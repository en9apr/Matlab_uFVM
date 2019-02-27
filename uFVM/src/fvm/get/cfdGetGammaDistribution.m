function Gamma = cfdGetGammaDistribution

theMesh = cfdGetMesh;
theNumberOfBElements =  theMesh.numberOfBElements;

Gamma_interior = ones(20,20);
Gamma_interior(4:10,11:17) = 1e2;
Gamma_interior(11:17,4:10) = 1e3;
Gamma_interior = Gamma_interior(:);

Gamma_boundary = ones(theNumberOfBElements,1);

Gamma = [Gamma_interior; Gamma_boundary];