function nbPatches = cfdGetNumberOfPatches
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theMesh = cfdGetMesh;
nbPatches = theMesh.numberOfPatches;