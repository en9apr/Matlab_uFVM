function phi_f = cfdInterpolateFromNodesToFaces(phiNodes)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmNodes = theMesh.nodes;

numberOfFaces = theMesh.numberOfFaces;

phi_f=zeros(numberOfFaces,1);

for iFace=1:numberOfFaces
    
    theFace = fvmFaces(iFace);
    iNodes = theFace.iNodes;
    
    C = theFace.centroid;
    
    localSumOfInverseDistance=0;
    localPhi=0;
    for iNode=iNodes
        theNode = fvmNodes(iNode);
        N=theNode.centroid;
        d=cfdMagnitude(C-N);
        
        localPhi = localPhi + phiNodes(iNode)/d;
        localSumOfInverseDistance = localSumOfInverseDistance+1/d;
    end
    localPhi = localPhi/localSumOfInverseDistance;
    
    phi_f(iFace) = localPhi;   
end