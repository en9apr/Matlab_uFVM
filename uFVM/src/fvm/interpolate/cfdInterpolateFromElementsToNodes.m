function phi_n = cfdInterpolateFromElementsToNodes(phi)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theMesh = cfdGetMesh;
fvmNodes = theMesh.nodes;
fvmElements = theMesh.elements;
fvmFaces = theMesh.faces;

numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

numberOfElements = theMesh.numberOfElements;
numberOfNodes = theMesh.numberOfNodes;
phi_n = zeros(numberOfNodes,1);

for iNode = 1:numberOfNodes
    theNode = fvmNodes(iNode);
    N = theNode.centroid;

    localPhiNode=0;
    localInverseDistanceSum = 0;

    if(isempty(theNode.iFaces(theNode.iFaces>numberOfInteriorFaces)))
        localElementIndices = theNode.iElements;

        for iElement = localElementIndices
            theElement = fvmElements(iElement);
            C = theElement.centroid;

            d = cfdMagnitude(N-C);
            localPhi = phi(iElement);

            localPhiNode = localPhiNode + localPhi/d;
            localInverseDistanceSum = localInverseDistanceSum + 1/d;       
        end
    else 
        localBFacesIndices = theNode.iFaces(theNode.iFaces>numberOfInteriorFaces);
        for iBFace = localBFacesIndices
            theFace = fvmFaces(iBFace);
            C = theFace.centroid;
            iBElement = numberOfElements+(iBFace-numberOfInteriorFaces);

            d = cfdMagnitude(N-C);
            localPhi = phi(iBElement);

            localPhiNode = localPhiNode + localPhi/d;
            localInverseDistanceSum = localInverseDistanceSum + 1/d;       
        end
        
    end
    localPhiNode = localPhiNode/localInverseDistanceSum;
    %
    %
    phi_n(iNode) = localPhiNode;
    %
end
