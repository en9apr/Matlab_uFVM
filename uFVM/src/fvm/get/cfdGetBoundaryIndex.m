function theBIndex = cfdGetBoundaryIndex(theBoundaryName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theBIndex=0;
theMesh= cfdGetMesh;
theNumberOfBoundaries = theMesh.numberOfBoundaries;
for iBoundary = 1:theNumberOfBoundaries
    if(strcmp(theMesh.boundaries(iBoundary).userName,theBoundaryName)==1)

        theBIndex= iBoundary;
        break;
    end
    
end
