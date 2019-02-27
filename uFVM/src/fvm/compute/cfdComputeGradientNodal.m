function phiGrad = cfdComputeGradientNodal(phi,theMesh)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%
if(nargin<2)
    theMesh = cfdGetMesh;
end

fvmElements = theMesh.elements;

%
%=====================================================
% INTERIOR Gradients 
%=====================================================
theNumberOfElements = theMesh.numberOfElements;
theNumberOfBElements = theMesh.numberOfBElements;
theNumberOfInteriorFaces = theMesh.numberOfInteriorFaces;
%
phiNodes = cfdInterpolateFromElementsToNodes(phi);
phi_f = cfdInterpolateFromNodesToFaces(phiNodes);
%
phiGrad = zeros(3,theNumberOfElements+theNumberOfBElements);
%-----------------------------------------------------
% INTERIOR FACES contribution to gradient 
%-----------------------------------------------------
fvmFaces = theMesh.faces;
% interpolate phi to faces
%
for iFace=1:theNumberOfInteriorFaces
   %
   theFace = fvmFaces(iFace);
   %
   iElement1 = theFace.iOwner;
   iElement2 = theFace.iNeighbour;
   %
   Sf = theFace.Sf;
   %
   %
   phiGrad(:,iElement1) = phiGrad(:,iElement1) + phi_f(iFace)*Sf;
   phiGrad(:,iElement2) = phiGrad(:,iElement2) - phi_f(iFace)*Sf;
   
end


%=====================================================
% BOUNDARY FACES contribution to gradient 
%=====================================================
for iBPatch=1:theNumberOfBElements
    %
    iBFace = theNumberOfInteriorFaces+iBPatch;
    iBElement = theNumberOfElements+iBPatch;
    theFace = fvmFaces(iBFace);
    %
    iElement1 = theFace.iOwner;
    %
    Sb = theFace.Sf;
    phi_b = phi(iBElement);
    %
    phiGrad(:,iElement1) = phiGrad(:,iElement1) + phi_b*Sb;
    
end


%-----------------------------------------------------
% Get Average Gradient by dividing with element volume 
%-----------------------------------------------------
for iElement =1:theNumberOfElements
   theElement = fvmElements(iElement);
   phiGrad(:,iElement) = phiGrad(:,iElement)/theElement.volume;
end


%-----------------------------------------------------
% Set boundary Gradient equal to associated element
% Gradient
%-----------------------------------------------------
for iBPatch = 1:theNumberOfBElements
    iBElement = iBPatch+theNumberOfElements;
    iBFace = iBPatch+theNumberOfInteriorFaces;
    theBFace = fvmFaces(iBFace);
    iOwner = theBFace.iOwner;
    phiGrad(:,iBElement) = phiGrad(:,iOwner);
end


