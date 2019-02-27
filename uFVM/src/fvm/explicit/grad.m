function phiGrad = grad(phi)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function computes the gradient of a scalar/vector field
%--------------------------------------------------------------------------

% Interior Gradients

% Get mesh info
theMesh = cfdGetMesh;
theNumberOfComponents = size(phi, 2);

iFaces = 1:theMesh.numberOfInteriorFaces;
iBFaces = theMesh.numberOfInteriorFaces+1:theMesh.numberOfFaces;

iElements = 1:theMesh.numberOfElements;
iBElements = theMesh.numberOfElements+1:theMesh.numberOfElements+theMesh.numberOfBFaces;

iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';

Sf = [theMesh.faces(iFaces).Sf]';
gf = [theMesh.faces(iFaces).gf]';

% Initialize
phiGrad = zeros(theMesh.numberOfElements+theMesh.numberOfBElements, 3, theNumberOfComponents);

% loop over interior faces and add contribution
for iComponent=1:theNumberOfComponents   
    phi_f = gf.*phi(iNeighbours,iComponent) + (1-gf).*phi(iOwners,iComponent);
    for iFace=iFaces
        phiGrad(iOwners(iFace),:,iComponent)     = phiGrad(iOwners(iFace),:,iComponent)     + phi_f(iFace)*Sf(iFace,:);
        phiGrad(iNeighbours(iFace),:,iComponent) = phiGrad(iNeighbours(iFace),:,iComponent) - phi_f(iFace)*Sf(iFace,:);
    end    
end

% Boundary faces contribution to gradient
iBOwners = [theMesh.faces(iBFaces).iOwner]';
phi_b = phi(iBElements,:);
Sb = [theMesh.faces(iBFaces).Sf]';
for iComponent=1:theNumberOfComponents
    for k=1:theMesh.numberOfBFaces
        phiGrad(iBOwners(k),:,iComponent) = phiGrad(iBOwners(k),:,iComponent) + phi_b(k)*Sb(k,:);
    end
end

% Get Average Gradient by dividing with element volume
volumes = [theMesh.elements(iElements).volume]';
for iComponent=1:theNumberOfComponents
    for iElement =1:theMesh.numberOfElements
        phiGrad(iElement,:,iComponent) = phiGrad(iElement,:,iComponent)/volumes(iElement);
    end
end

% Set boundary Gradient equal to associated element gradient
phiGrad(iBElements,:,:) = phiGrad(iBOwners,:,:);


% Boundary Gradients
theNumberOfPatches = theMesh.numberOfBoundaries;

for iPatch=1:theNumberOfPatches
    % find the Physical Type
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    % Wall
    if(strcmp(thePhysicalType,'wall'))
        phiGrad = cfdUpdateWallGradients(iPatch,phi,phiGrad);
    end
    % Inlet
    if(strcmp(thePhysicalType,'inlet'))
        phiGrad = cfdUpdateInletGradients(iPatch,phi,phiGrad);
    end
    % Outlet
    if(strcmp(thePhysicalType,'outlet'))
        phiGrad = cfdUpdateOutletGradients(iPatch,phi,phiGrad);
    end
    % Symmetry
    if(strcmp(thePhysicalType,'symmetry'))
        phiGrad = cfdUpdateSymmetryGradients(iPatch,phi,phiGrad);
    end
    
    % Empty
    if(strcmp(thePhysicalType,'empty'))
        phiGrad = cfdUpdateSymmetryGradients(iPatch,phi,phiGrad);
    end
    
    % Patch
    if(strcmp(thePhysicalType,'patch'))
        phiGrad = cfdUpdateWallGradients(iPatch,phi,phiGrad);
    end
        
end

end



%=====================================================
% WALL Gradients
%=====================================================

function grad = cfdUpdateWallGradients(iPatch, phi, grad)
%
%
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmElements = theMesh.elements;
theBoundary = theMesh.boundaries(iPatch);
%
theNumberOfComponents = 1;
if size(phi, 2)>1
    theNumberOfComponents = 3;
end

numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
endBFace = startBFace+numberOfBFaces-1;

startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


% patch gradient
grad_b = zeros(numberOfBFaces,3,theNumberOfComponents);
%
for iComponent = 1:theNumberOfComponents
    
    % loop over patchFaces
    for iBPatch=1:numberOfBFaces
        iBFace = startBFace-1+iBPatch;
        iBElement = startBElement-1+iBPatch;
        theFace = fvmFaces(iBFace);
        iOwner = theFace.iOwner;
        theElement = fvmElements(iOwner);
        d = (theFace.centroid - theElement.centroid)';
        dmag = norm(d);
        e = d/dmag;
        grad_b(iBPatch,:,iComponent) = grad(iOwner,:,iComponent) - (grad(iOwner,:,iComponent)*e')*e + ((phi(iBElement,iComponent) - phi(iOwner,iComponent))/dmag)*e;
    end
end
%
grad(startBElement:endBElement,:,:) = grad_b;
%
end


%=====================================================
% INLET Gradients
%=====================================================


function grad = cfdUpdateInletGradients(iPatch, phi, grad)
%
%
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmElements = theMesh.elements;
theBoundary = theMesh.boundaries(iPatch);
%
theNumberOfComponents = 1;
if size(phi, 2)>1
    theNumberOfComponents = 3;
end
%

numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
endBFace = startBFace+numberOfBFaces-1;

startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


% patch gradient
grad_b = zeros(numberOfBFaces,3,theNumberOfComponents);

for iComponent = 1:theNumberOfComponents
    % loop over patchFaces
    for iBPatch=1:numberOfBFaces
        iBFace = startBFace-1+iBPatch;
        iBElement = startBElement-1+iBPatch;
        theFace = fvmFaces(iBFace);
        iOwner = theFace.iOwner;
        theElement = fvmElements(iOwner);
        d = (theFace.centroid - theElement.centroid)';
        dmag = norm(d);
        e = d/dmag;
        
        grad_b(iBPatch,:,iComponent) = grad(iOwner,:,iComponent) - (grad(iOwner,:,iComponent)*e')*e + ((phi(iBElement) - phi(iOwner))/dmag)*e;
    end
end
%
grad(startBElement:endBElement,:,:) = grad_b;
end
%


%=====================================================
% OUTLET Gradients
%=====================================================
function grad = cfdUpdateOutletGradients(iPatch, phi, grad)
%
%
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmElements = theMesh.elements;
theBoundary = theMesh.boundaries(iPatch);
%
%
theNumberOfComponents = 1;
if size(phi, 2)>1
    theNumberOfComponents = 3;
end
%

numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
theEndBFace = startBFace+numberOfBFaces-1;

startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


% patch gradient
grad_b = zeros(numberOfBFaces,3,theNumberOfComponents);
for iComponent = 1:theNumberOfComponents
    % loop over patchFaces
    for iBPatch=1:numberOfBFaces
        iBFace = startBFace-1+iBPatch;
        iBElement = startBElement-1+iBPatch;
        theFace = fvmFaces(iBFace);
        iOwner = theFace.iOwner;
        theElement = fvmElements(iOwner);
        d = (theFace.centroid - theElement.centroid)';
        dmag = norm(d);
        e = d/dmag;
        
        grad_b(iBPatch,:,iComponent) = grad(iOwner,:,iComponent) - (grad(iOwner,:,iComponent)*e')*e + ((phi(iBElement) - phi(iOwner))/dmag)*e;
    end
end
%
grad(startBElement:endBElement,:,:) = grad_b;
%
end

%=====================================================
% SYMMETRY Gradients
%=====================================================
function grad = cfdUpdateSymmetryGradients(iPatch,phi, grad)

end
