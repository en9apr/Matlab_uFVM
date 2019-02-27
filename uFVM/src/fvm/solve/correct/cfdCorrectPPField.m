function cfdCorrectPPField
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%
%=============================
% Correct PP @ Interior
%=============================
cfdCorrectPPInterior;
%=============================
% Correct PP @ Patches
%=============================
cfdCorrectPPBoundaryPatches;
%=============================
% Updage PP Gradient
%=============================
cfdUpdateGradient('PP');
%
%
end

%===================================================
% Correct Interior
%===================================================
function cfdCorrectPPInterior
%
thePPField = cfdGetMeshField('PP');
%
iFixedPressureIndex = cfdGetFixedElement;
%
pp  = thePPField.phi;
dphi = cfdGetDPhi;
%=======================================================================
% Loop over all cells and update the field: phi@n+1 = phi@n + dphi@n+1
%=======================================================================
theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;
iElements = 1:theNumberOfElements;
%
pp(iElements) = dphi;
dphi = zeros(theNumberOfElements,1);

ppc = 0;
if(iFixedPressureIndex>0)
    ppc = pp(iFixedPressureIndex);
end
cfdSetDPhi(dphi);
%
thePPField.phi(iElements) = pp(iElements) - ppc;
cfdSetMeshField(thePPField);
end

%===================================================
% Correct Boundary Patches
%===================================================
function cfdCorrectPPBoundaryPatches

theMesh = cfdGetMesh;

thePressureEquation = cfdGetModel('p');
%
theNumberOfPatches = theMesh.numberOfPatches;
for iPatch=1:theNumberOfPatches
    
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    theBCType = thePressureEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if strcmp(thePhysicalType,'wall')
        if strcmp(theBCType,'slip') || strcmp(theBCType,'noSlip') || strcmp(theBCType,'zeroGradient')
            cfdCorrectPPWallBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'inlet') || strcmp(theBCType,'zeroGradient')
            cfdCorrectPPInletInletBC(iPatch);
        elseif(strcmp(theBCType,'fixedValue'))
            cfdCorrectPPInletFixedValueBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalType,'outlet')
        if strcmp(theBCType,'outlet')
            cfdCorrectPPOutletOutletBC(iPatch);
        elseif strcmp(theBCType,'fixedValue')
            cfdCorrectPPOutletFixedValueBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % SYMMETRY
        %
    elseif strcmp(thePhysicalType,'symmetry')
        if strcmp(theBCType,'symmetry')
            cfdCorrectPPSymmetryBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % EMPTY
        %
    elseif strcmp(thePhysicalType,'empty')
        if strcmp(theBCType,'empty')
            cfdCorrectPPEmptyBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
    else
        error('Physical Boundary condition Not Implemented');
        
    end
    %
end
%
end


%===================================================
% WALL-noslip Wall
%===================================================
function cfdCorrectPPWallBC(iPatch)
%
thePPField = cfdGetMeshField('PP');
phi = thePPField.phi;

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theBoundary = theMesh.boundaries(iPatch);

startBFace = theBoundary.startFace;
endBFace = startBFace+theBoundary.numberOfBFaces-1;
iBFaces = startBFace:endBFace;

startBElement = theMesh.numberOfElements+startBFace-numberOfInteriorFaces;
endBElement = startBElement+theBoundary.numberOfBFaces-1;
iBElements = startBElement:endBElement;
%
iOwners = [theMesh.faces(iBFaces).iOwner];
phi(iBElements) = phi(iOwners);

%
thePPField.phi = phi;
cfdSetMeshField(thePPField);
end

%===================================================
% INLET-INLET
%===================================================
function cfdCorrectPPInletInletBC(iPatch)
%
thePPField = cfdGetMeshField('PP');
phi = thePPField.phi;

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theBoundary = theMesh.boundaries(iPatch);

startBFace = theBoundary.startFace;
endBFace = startBFace+theBoundary.numberOfBFaces-1;
iBFaces = startBFace:endBFace;

startBElement = theMesh.numberOfElements+startBFace-numberOfInteriorFaces;
endBElement = startBElement+theBoundary.numberOfBFaces-1;
iBElements = startBElement:endBElement;
%
iOwners = [theMesh.faces(iBFaces).iOwner];
phi(iBElements) = phi(iOwners);

%
thePPField.phi = phi;
cfdSetMeshField(thePPField);
end


%===================================================
% INLET-specifiedValue
%===================================================
function cfdCorrectPPInletFixedValueBC(iPatch)
%
thePPField = cfdGetMeshField('PP');
phi = thePPField.phi;

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theBoundary = theMesh.boundaries(iPatch);

startBElement = theMesh.numberOfElements+startBFace-numberOfInteriorFaces;
endBElement = startBElement+theBoundary.numberOfBFaces-1;
iBElements = startBElement:endBElement;
%
phi(iBElements) = 0;

%
thePPField.phi = phi;
cfdSetMeshField(thePPField);
end


%===================================================
% OUTLET-Outlet
%===================================================
function cfdCorrectPPOutletOutletBC(iPatch)
%

thePPField = cfdGetMeshField('PP');
phi = thePPField.phi;

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theBoundary = theMesh.boundaries(iPatch);

startBFace = theBoundary.startFace;
endBFace = startBFace+theBoundary.numberOfBFaces-1;
iBFaces = startBFace:endBFace;

startBElement = theMesh.numberOfElements+startBFace-numberOfInteriorFaces;
endBElement = startBElement+theBoundary.numberOfBFaces-1;
iBElements = startBElement:endBElement;
%
iOwners = [theMesh.faces(iBFaces).iOwner];
phi(iBElements) = phi(iOwners);

%
thePPField.phi = phi;
cfdSetMeshField(thePPField);
end

%===================================================
% OUTLET-specifiedValue
%===================================================
function cfdCorrectPPOutletFixedValueBC(iPatch)
%

thePPField = cfdGetMeshField('PP');
phi = thePPField.phi;

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theBoundary = theMesh.boundaries(iPatch);

startBFace = theBoundary.startFace;

startBElement = theMesh.numberOfElements+startBFace-numberOfInteriorFaces;
endBElement = startBElement+theBoundary.numberOfBFaces-1;
iBElements = startBElement:endBElement;
%
phi(iBElements) = 0;

%
thePPField.phi = phi;
cfdSetMeshField(thePPField);
end


%===================================================
% SYMMETRY
%===================================================
function cfdCorrectPPSymmetryBC(iPatch)
%

thePPField = cfdGetMeshField('PP');
phi = thePPField.phi;

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theBoundary = theMesh.boundaries(iPatch);

startBFace = theBoundary.startFace;
endBFace = startBFace+theBoundary.numberOfBFaces-1;
iBFaces = startBFace:endBFace;

startBElement = theMesh.numberOfElements+startBFace-numberOfInteriorFaces;
endBElement = startBElement+theBoundary.numberOfBFaces-1;
iBElements = startBElement:endBElement;
%
iOwners = [theMesh.faces(iBFaces).iOwner];
phi(iBElements) = phi(iOwners);

%
thePPField.phi = phi;
cfdSetMeshField(thePPField);
end

%===================================================
% Empty
%===================================================
function cfdCorrectPPEmptyBC(iPatch)
%

thePPField = cfdGetMeshField('PP');
phi = thePPField.phi;

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theBoundary = theMesh.boundaries(iPatch);

startBFace = theBoundary.startFace;
endBFace = startBFace+theBoundary.numberOfBFaces-1;
iBFaces = startBFace:endBFace;

startBElement = theMesh.numberOfElements+startBFace-numberOfInteriorFaces;
endBElement = startBElement+theBoundary.numberOfBFaces-1;
iBElements = startBElement:endBElement;
%
iOwners = [theMesh.faces(iBFaces).iOwner];
phi(iBElements) = phi(iOwners);

%
thePPField.phi = phi;
cfdSetMeshField(thePPField);
end



