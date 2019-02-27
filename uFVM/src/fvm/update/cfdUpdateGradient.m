function cfdUpdateGradient(theFieldName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function updates gradient of the field
%--------------------------------------------------------------------------
%
theField = cfdGetModel(theFieldName);

gradientType = 'Gauss linear';
if(isfield(theField,'gradientType'))
    gradientType = theField.gradientType;
end

if(strcmp(gradientType,'Gauss linear'))
    cfdUpdateGradientGauss0(theFieldName);
    %
elseif(strcmp(gradientType,'NODAL') || strcmp(gradientType,'Nodal'))
    cfdUpdateGradientNodal(theFieldName);
    %
elseif(strcmp(gradientType,'LEAST SQUARE') || strcmp(gradientType,'LEASTSQUARE') || strcmp(gradientType,'LS0') || strcmp(gradientType,'LS'))
    cfdUpdateGradientLS0(theFieldName);
    disp('****LS***')
    %
elseif(strcmp(gradientType,'LEAST SQUARE-2') || strcmp(gradientType,'LEASTSQUARE-2'))
    cfdUpdateGradientLeastSquare(theFieldName);
    %
else
    cfdUpdateGradientGauss0(theFieldName);
end


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
        cfdUpdateWallGradients(iPatch,theFieldName);
    end
    % Inlet
    if(strcmp(thePhysicalType,'inlet'))
        cfdUpdateInletGradients(iPatch,theFieldName);
    end
    % Outlet
    if(strcmp(thePhysicalType,'outlet'))
        cfdUpdateOutletGradients(iPatch,theFieldName);
    end
    % Symmetry
    if(strcmp(thePhysicalType,'symmetry'))
        cfdUpdateSymmetryGradients(iPatch,theFieldName);
    end
    % Patch
    if(strcmp(thePhysicalType,'patch'))
        cfdUpdateWallGradients(iPatch,theFieldName);
    end
    
    
end

end



%=====================================================
% WALL Gradients
%=====================================================


function cfdUpdateWallGradients(iPatch,theFieldName)
%
%
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmElements = theMesh.elements;
theBoundary = theMesh.boundaries(iPatch);
%
theMeshField = cfdGetMeshField(theFieldName);
theType = theMeshField.type;

theNumberOfComponents = 1;
if(strcmp(theType,'Vector'))
    theNumberOfComponents = 3;
end
%
phi = theMeshField.phi;
%

numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
endBFace = startBFace+numberOfBFaces-1;

startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


% patch gradient
grad = theMeshField.phiGradient;
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
theMeshField.phiGradient(startBElement:endBElement,:,:) = grad_b;
cfdSetMeshField(theMeshField);
%
end


%=====================================================
% INLET Gradients
%=====================================================


function cfdUpdateInletGradients(iPatch,theFieldName)
%
%
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmElements = theMesh.elements;
theBoundary = theMesh.boundaries(iPatch);
%
%
theMeshField = cfdGetMeshField(theFieldName);
theType = theMeshField.type;
%
theNumberOfComponents = 1;
if(strcmp(theType,'Vector'))
    theNumberOfComponents= 3;
end
%
phi = theMeshField.phi;
%

numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
endBFace = startBFace+numberOfBFaces-1;

startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


% patch gradient
grad = theMeshField.phiGradient;
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
theMeshField.phiGradient(startBElement:endBElement,:,:) = grad_b;
cfdSetMeshField(theMeshField);
end
%


%=====================================================
% OUTLET Gradients
%=====================================================
function cfdUpdateOutletGradients(iPatch,theFieldName)
%
%
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
fvmElements = theMesh.elements;
theBoundary = theMesh.boundaries(iPatch);
%
%
theMeshField = cfdGetMeshField(theFieldName);
theType = theMeshField.type;
%
theNumberOfComponents = 1;
if(strcmp(theType,'Vector'))
    theNumberOfComponents= 3;
end
%
phi = theMeshField.phi;
%

numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
theEndBFace = startBFace+numberOfBFaces-1;

startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


% patch gradient
grad = theMeshField.phiGradient;
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
theMeshField.phiGradient(startBElement:endBElement,:,:) = grad_b;
cfdSetMeshField(theMeshField);
%
end

%=====================================================
% SYMMETRY Gradients
%=====================================================
function cfdUpdateSymmetryGradients(iPatch,theFieldName)

end
