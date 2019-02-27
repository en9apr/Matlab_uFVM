function theIndices = cfdGetFaceIndicesForBoundaryIndex(theBoundaryIndex)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theBoundary = cfdGetBoundary(theBoundaryIndex);

theNumberOfBFaces = theBoundary.numberOfBFaces;
theStartFace = theBoundary.startFace;

theIndices = [theStartFace:theStartFace+theNumberOfBFaces-1];

end