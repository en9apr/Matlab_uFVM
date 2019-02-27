function theBPatchUserNames = cfdGetBPatchUserNames
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theMesh = cfdGetMesh;

theBPatchUserNames = {};
theNumberOfBPatches = theMesh.numberOfBoundaries;
for iBPatch=1:theNumberOfBPatches
    theBoundary = theMesh.boundaries(iBPatch);
    theBPatchUserNames{iBPatch} = theBoundary.userName;
    
end