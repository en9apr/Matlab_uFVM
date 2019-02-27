function cfdSetupPhysicalConditions(thePCCells);

global CFDEnv;

theMesh = cfdGetMesh;
theMesh.boundaryTypes = fvmDefineBoundaryTypes(thePCCells);

CFDEnv.mesh = theMesh;