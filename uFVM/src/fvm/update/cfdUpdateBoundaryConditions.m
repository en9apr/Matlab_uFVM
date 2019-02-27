function theField = cfdUpdateBoundaryConditions(theFieldName)

%
%===================================================
%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%

theMesh = cfdGetMesh;
theField = cfdGetModel(theFieldName);
theNumberOfPatches = theMesh.numberOfBoundaries;

for iPatch=1:theNumberOfPatches
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    theBCType = theField.bcs{iPatch}.type;
    %
    % WALL
    %
    if(strcmp(thePhysicalType,'wall'))
        if(strcmp(theBCType,'fixedValue'))
            cfdUpdateWallFixedValueBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'fixedGradient'))
            cfdUpdateWallFixedGradientBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'zeroGradient'))
            cfdUpdateWallZeroGradientBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'hybrid') || strcmp(theBCType,'Hybrid'))
            cfdUpdateWallHybridBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'noSlip'))
            cfdUpdateWallNoslipBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'slip'))
            cfdUpdateWallslipBC(iPatch,theFieldName);
        else
            error([theBCType 'Wall bc not defined']);
        end
        %
        % PATCH
        %
    elseif(strcmp(thePhysicalType,'patch'))
        if(strcmp(theBCType,'fixedValue'))
            cfdUpdatePatchFixedValueBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'fixedGradient'))
            cfdUpdateOutletFixedValueBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'zeroFlux')||strcmp(theBCType,'zeroGradient'))
            cfdUpdateOutletFixedValueBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'outlet'))
            cfdUpdateOutletZeroGradientBC(iPatch,theFieldName);
        else
            error([theBCType 'Wall bc not defined']);
        end
        %
        % INLET
        %
    elseif(strcmp(thePhysicalType,'inlet'))
        if(strcmp(theBCType,'fixedValue'))
            cfdUpdateInletFixedValueBC(iPatch,theFieldName);
        elseif(strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'Inlet'))
            cfdUpdateInletZeroGradientBC(iPatch,theFieldName);
        else
            error('Inlet bc not defined');
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalType,'outlet')
        if strcmp(theBCType,'fixedValue')
            cfdUpdateOutletFixedValueBC(iPatch,theFieldName);            
        elseif strcmp(theBCType,'outlet') || strcmp(theBCType,'zeroGradient')
            cfdUpdateOutletZeroGradientBC(iPatch,theFieldName);            
        else
            error([theBCType 'Outlet bc not defined']);
        end
        %
        % SYMMETRY
        %
    elseif(strcmp(thePhysicalType,'symmetry'))
        if(strcmp(theBCType,'symmetry') || strcmp(theBCType,'Symmetry'))
            cfdUpdateSymmetrySymmetryBC(iPatch,theFieldName);
        else
            error('Symmetry bc not defined');
        end
        
    elseif(strcmp(thePhysicalType,'empty')|| strcmp(theBCType,'Empty'))
        if(strcmp(theBCType,'empty') || strcmp(theBCType,'Empty'))
            cfdUpdateEmptyBC(iPatch,theFieldName);
        else
            error('Empty bc not defined');
        end
    else
        error([thePhysicalType '<<<< Physical Condition bc not defined']);
        
    end
                
end

end



%===================================================
% WALL-fixedValue
%===================================================
function cfdUpdateWallFixedValueBC(iPatch,theFieldName)
%
theMesh = cfdGetMesh;
theField = cfdGetModel(theFieldName);

theBCValue = theField.bcs{iPatch}.value;
% boundary phi values
% update the boundary values

theType = theField.type;
theLocale = cfdGetLocaleForPatch(iPatch);

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


if strcmp(theType,'Scalar')
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement) = phi_b .* ones(numberOfBFaces,1);
    cfdSetMeshField(theMeshField);
    
elseif strcmp(theType,'Vector')
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement,:) = phi_b .* ones(numberOfBFaces,3);
    cfdSetMeshField(theMeshField);
    
end

end

%===================================================
% WALL-fixedGradient
%===================================================
function cfdUpdateWallFixedGradientBC(iPatch,theFieldName)

% Get phi
theMeshField = cfdGetMeshField(theFieldName);
phi = theMeshField.phi;

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iBElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iBElementEnd = iBElementStart+numberOfBFaces-1;
iBElements = iBElementStart:iBElementEnd;

% Compute phiFlux
theEquation = cfdGetModel(theFieldName);
formula = theEquation.bcs{iPatch}.value;
theLocale = ['BPatch' num2str(iPatch)];

% Get gamma
theDiffusionTerm = cfdGetTermInEquation(theFieldName,'Diffusion');
theGammaField = cfdGetMeshField(theDiffusionTerm.coefficientName);
gamma = theGammaField.phi(iBElements);

% Calculate heat flux
normGrad_b = cfdComputeFormulaAtLocale(formula, theLocale);
phiFlux = - gamma .* normGrad_b;

% Get other info
gDiff = [theMesh.faces(iBFaces).gDiff]';
area = [theMesh.faces(iBFaces).area]';
%
iOwners = [theMesh.faces(iBFaces).iOwner];

% Update boundary phi
phi(iBElements) = (-phiFlux.*area - gamma .* gDiff .* phi(iOwners)) ./ (gamma .* gDiff);

% Update meshfield
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end

%===================================================
% WALL-zeroFlux
%===================================================
function cfdUpdateWallZeroGradientBC(iPatch,theFieldName)

end
%===================================================
% WALL-Hybrid
%===================================================
function cfdUpdateWallHybridBC(iPatch,theFieldName)
%

end


%===================================================
% WALL-noSlip
%===================================================
function cfdUpdateWallNoslipBC(iPatch,theFieldName)
%
theMesh=cfdGetMesh;
theField = cfdGetModel(theFieldName);
theBCValue = theField.bcs{iPatch}.value;
% boundary phi values
% update the boundary values
theType = theField.type;
theLocale = cfdGetLocaleForPatch(iPatch);

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;

theFieldBaseName = cfdGetBaseName(theField.name);

if strcmp(theField.name,'p')
    
elseif strcmp(theFieldBaseName,'TKE')
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    
    tke = theField.phi;
    tke_b = theField.phiPatches{iPatch};
    
    theNumberOfPatchElements = theMesh.patchSizes(iPatch);
    patchFacesIndices = find(fvmCellsToArray({theMesh.faces(:).patchIndex})==iPatch);
    for k = theNumberOfPatchElements
        iBFace = patchFacesIndices(k);
        theFace = theMesh.faces(iBFace);
        iElement1= theFace.element1;
        tke_b(k) = tke(iElement1);
    end
    theField.phiPatches{iPatch} = tke_b;
    
elseif strcmp(theFieldBaseName,'TED')
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    
    
    ted = theField.phi;
    ted_b = theField.phiPatches{iPatch};
    
    theNumberOfPatchElements = theMesh.patchSizes(iPatch);
    patchFacesIndices = find(fvmCellsToArray({theMesh.faces(:).patchIndex})==iPatch);
    for k = theNumberOfPatchElements
        iBFace = patchFacesIndices(k);
        theFace = theMesh.faces(iBFace);
        iElement1= theFace.element1;
        ted_b(k) = ted(iElement1);
    end
    theField.phiPatches{iPatch} = ted_b;
    
elseif strcmp(theFieldBaseName,'U')
    
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement,:) = phi_b .* ones(numberOfBFaces,3);
    cfdSetMeshField(theMeshField);
end

end


%===================================================
% WALL-slip
%===================================================
function cfdUpdateWallslipBC(iPatch,theFieldName)

end

%===================================================
% INLET-fixedValue
%===================================================
function cfdUpdateInletFixedValueBC(iPatch,theFieldName)
%
theMesh = cfdGetMesh;
theField = cfdGetModel(theFieldName);
theBCValue = theField.bcs{iPatch}.value;
% boundary phi values
% update the boundary values
theType = theField.type;
theLocale = cfdGetLocaleForPatch(iPatch);

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


if strcmp(theType,'Scalar')
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement) = phi_b .* ones(numberOfBFaces,1);
    cfdSetMeshField(theMeshField);    
elseif strcmp(theType,'Vector')    
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement,:) = phi_b .* ones(numberOfBFaces,3);
    cfdSetMeshField(theMeshField);
end


end

%===================================================
% INLET-Inlet
%===================================================
function cfdUpdateInletZeroGradientBC(iPatch,theFieldName)

end
%===================================================
% OUTLET-fixedValue
%===================================================
function cfdUpdateOutletFixedValueBC(iPatch,theFieldName)

%
theMesh=cfdGetMesh;
theField = cfdGetModel(theFieldName);
theBCValue = theField.bcs{iPatch}.value;
% boundary phi values
% update the boundary values

theType = theField.type;
theLocale = cfdGetLocaleForPatch(iPatch);

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;


if strcmp(theType,'Scalar')
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement) = phi_b .* ones(numberOfBFaces,1);
    cfdSetMeshField(theMeshField);    
elseif strcmp(theType,'Vector')    
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement,:) = phi_b .* ones(numberOfBFaces,3);
    cfdSetMeshField(theMeshField);
end


end
%===================================================
% OUTLET-Outlet
%===================================================
function cfdUpdateOutletZeroGradientBC(iPatch,theFieldName)

end

%===================================================
% SYMMETRY-Symmetry
%===================================================
function cfdUpdateSymmetrySymmetryBC(iPatch,theFieldName)

end

%===================================================
% EMPTY-empty
%===================================================
function cfdUpdateEmptyBC(iPatch,theFieldName)

end

%===================================================
% PATCH-fixedValue
%===================================================
function cfdUpdatePatchFixedValueBC(iPatch,theFieldName)
%

theMesh = cfdGetMesh;
theField = cfdGetModel(theFieldName);
theBCValue = theField.bcs{iPatch}.value;

theType = theField.type;
theLocale = cfdGetLocaleForPatch(iPatch);

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

startBFace = theBoundary.startFace;
startBElement = startBFace-numberOfInteriorFaces+numberOfElements;
endBElement = startBElement+numberOfBFaces-1;

if strcmp(theType,'Scalar')
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement) = phi_b .* ones(numberOfBFaces,1);
    cfdSetMeshField(theMeshField);    
elseif strcmp(theType,'Vector')   
    phi_b = cfdComputeFormulaAtLocale(theBCValue,theLocale,theType);
    theMeshField = cfdGetMeshField(theField.name);
    theMeshField.phi(startBElement:endBElement,:) = phi_b .* ones(numberOfBFaces,3);
    cfdSetMeshField(theMeshField);
end


end
