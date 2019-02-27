function theBoundary = cfdGetBoundary(thePatchIndex)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(thePatchIndex);

end