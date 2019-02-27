function theIndices = cfdGetElementIndicesForBoundaryIndex(theBoundaryIndex)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

thePreviousPatchesNumberOfBElements = 0;
for iPatch=1:theBoundaryIndex-1
    theBoundary = cfdGetBoundary(iPatch);
    thePreviousPatchesNumberOfBElements = thePreviousPatchesNumberOfBElements+theBoundary.numberOfBFaces;
end

theBoundary = cfdGetBoundary(theBoundaryIndex);
theNumberOfBElements = theBoundary.numberOfBFaces;

theNumberOfElements = cfdGetNumberOfElements;
theStartElement = theNumberOfElements+thePreviousPatchesNumberOfBElements+1;

theIndices = [theStartElement:theStartElement+theNumberOfBElements-1];

end