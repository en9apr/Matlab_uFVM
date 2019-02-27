function cfdCorrectPressureEquation
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%=============================
% Correct Pressure @ interior
%=============================
cfdCorrectPressureInterior;
%=============================
% Correct Pressure @ BC
%=============================
cfdCorrectPressureBoundaryPatches;
%=============================

end


%===================================================
% Correct Interior
%===================================================
function cfdCorrectPressureInterior

% Get mesh info
theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;
iElements = 1:theNumberOfElements;

% Get pressure field
thePressureEquation = cfdGetModel('p');
thePressureField = cfdGetMeshField('p');
p = thePressureField.phi;

% Get pressure correction field
thePPField  = cfdGetMeshField('PP');
pp = thePPField.phi;

% Get pressure under-relaxation
URFP = thePressureEquation.urf;
%
%
ppFixed = 0;
iFixedElement = cfdGetFixedElement;
if(iFixedElement>0)
    ppFixed = pp(iFixedElement);
end
pp = pp - ppFixed;

% Update pressure field with explicit under-relaxation
p(iElements) =  p(iElements) + URFP*pp(iElements);

thePressureField.phi = p;
cfdSetMeshField(thePressureField);
end


%===================================================
% Correct Boundary Patches
%===================================================
function cfdCorrectPressureBoundaryPatches

theMesh = cfdGetMesh;
theNumberOfPatches = theMesh.numberOfPatches;
thePressureEquation = cfdGetModel('p');

for iPatch=1:theNumberOfPatches
    
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    theBCType = thePressureEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if strcmp(thePhysicalType,'wall')
        if strcmp(theBCType,'slip') || strcmp(theBCType,'noSlip') || strcmp(theBCType,'zeroGradient')
            cfdCorrectPressureWallBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'inlet') || strcmp(theBCType,'zeroGradient')
            cfdCorrectPressureInletInletBC(iPatch);
        elseif strcmp(theBCType,'fixedValue')
            cfdCorrectPressureInletFixedValueBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalType,'outlet')
        if strcmp(theBCType,'outlet')
            cfdCorrectPressureOutletOutletBC(iPatch);
        elseif(strcmp(theBCType,'fixedValue'))
            cfdCorrectPressureOutletFixedValueBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % SYMMETRY
        %
    elseif strcmp(thePhysicalType,'symmetry')
        if strcmp(theBCType,'symmetry')
            cfdCorrectPressureSymmetryBC(iPatch);
        else
            error([theBCType 'BC Condition not Implemented']);
        end
        %
        % EMPTY
        %
    elseif strcmp(thePhysicalType,'empty')
        cfdCorrectPressureEmptyBC(iPatch);
    else
        %
        % ERROR
        %
        error('Physical Condition in Pressure Correction not Implemented');
    end
    
end


end

%===================================================
% WALL
%===================================================
function cfdCorrectPressureWallBC(iPatch)
%
thePressureField = cfdGetMeshField('p');
p = thePressureField.phi;
% PatchFACES
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
p(iBElements) = p(iOwners);
%
thePressureField.phi = p;
cfdSetMeshField(thePressureField);
end

%===================================================
% INLET-INLET
%===================================================
function cfdCorrectPressureInletInletBC(iPatch)
%
thePressureField = cfdGetMeshField('p');
p = thePressureField.phi;
% PatchFACES
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
p(iBElements) = p(iOwners);
%
thePressureField.phi = p;
cfdSetMeshField(thePressureField);
end


%===================================================
% INLET-specifiedValue
%===================================================
function cfdCorrectPressureInletFixedValueBC(iPatch)
%


end


%===================================================
% OUTLET-Outlet
%===================================================
function cfdCorrectPressureOutletOutletBC(iPatch)
%
%
thePressureField = cfdGetMeshField('p');
p = thePressureField.phi;
% PatchFACES
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
p(iBElements) = p(iOwners);
%
thePressureField.phi = p;
cfdSetMeshField(thePressureField);
end

%===================================================
% OUTLET-specifiedValue
%===================================================
function cfdCorrectPressureOutletFixedValueBC(iPatch)


end
%===================================================
% Symmetry BC
%===================================================

function cfdCorrectPressureSymmetryBC(iPatch)

thePressureField = cfdGetMeshField('p');
p = thePressureField.phi;

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
p(iBElements) = p(iOwners);
%
thePressureField.phi = p;
cfdSetMeshField(thePressureField);
end

%===================================================
% Empty BC
%===================================================

function cfdCorrectPressureEmptyBC(iPatch)

thePressureField = cfdGetMeshField('p');
p = thePressureField.phi;

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
p(iBElements) = p(iOwners);
%
thePressureField.phi = p;
cfdSetMeshField(thePressureField);
end
