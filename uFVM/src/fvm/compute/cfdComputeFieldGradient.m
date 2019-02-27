function phiGrad = cfdComputeFieldGradient(phi)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function updates gradient of the field
%--------------------------------------------------------------------------

% Compute green gauss gradient
phiGrad = cfdComputeGradientGauss0(phi');

%=====================================================
% BOUNDARY Gradients
%=====================================================
theMesh = cfdGetMesh;

theNumberOfPatches = theMesh.numberOfBoundaries;

for iPatch=1:theNumberOfPatches
    % find the Physical Type
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    % Wall
    if(strcmp(thePhysicalType,'wall'))
        % find type of BC
        phiGrad = cfdUpdateWallGradients(iPatch,phi',phiGrad);
    end
    % Inlet
    if(strcmp(thePhysicalType,'inlet'))
        phiGrad = cfdUpdateInletGradients(iPatch,phi',phiGrad);
    end
    % Outlet
    if(strcmp(thePhysicalType,'outlet'))
        phiGrad = cfdUpdateOutletGradients(iPatch,phi',phiGrad);
    end
    % Symmetry
    if(strcmp(thePhysicalType,'symmetry'))
        phiGrad = cfdUpdateSymmetryGradients(iPatch,phi',phiGrad);
    end
    % Patch
    if(strcmp(thePhysicalType,'patch'))
        phiGrad = cfdUpdateWallGradients(iPatch,phi',phiGrad);
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
