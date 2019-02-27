function iFaces = cfdGetiFacesForPatch(iPatch)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


m=cfdGetMesh;

theBoundary = m.boundaries(iPatch);
iFaces = theBoundary.startFace:theBoundary.startFace+theBoundary.numberOfBFaces-1;
