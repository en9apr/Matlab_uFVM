function cfdCorrectScalarEquation(theEquationName,iComponent)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

if(nargin==1)
    iComponent = 1;
end

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;
theScalarEquation = cfdGetModel(theEquationName);
theScalarMeshField = cfdGetMeshField(theEquationName);

dphi = cfdGetDPhi;
%---------------------------------
% Correct Interior Field
%---------------------------------
theScalarMeshField.phi(1:numberOfElements,iComponent) = theScalarMeshField.phi(1:numberOfElements,iComponent) + dphi;
cfdSetMeshField(theScalarMeshField);

dphi = zeros(numberOfElements,1);
cfdSetDPhi(dphi);
%
%---------------------------------
% Correct Boundary patches
%---------------------------------
theNumberOfBoundaries = theMesh.numberOfBoundaries;
for iPatch=1:theNumberOfBoundaries
    % get Physical and Equation Boundary Conditions
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    theBCType = theScalarEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if strcmp(thePhysicalType,'wall')
        if strcmp(theBCType,'fixedValue')
            cfdCorrectWallFixedValueBC(iPatch,theEquationName,iComponent);
        elseif strcmp(theBCType,'fixedGradient')
            cfdCorrectWallFixedGradientBC(iPatch,theEquationName,iComponent);
        elseif strcmp(theBCType,'zeroGradient')
            cfdCorrectWallZeroGradientBC(iPatch,theEquationName,iComponent);
        elseif strcmp(theBCType,'Hybrid')
            cfdCorrectWallHybridBC(iPatch,theEquationName,iComponent);
        elseif strcmp(theBCType,'noSlip')
            cfdCorrectWallNoslipBC(iPatch,theEquationName,iComponent);
        elseif strcmp(theBCType,'slip')
            cfdCorrectWallSlipBC(iPatch,theEquationName,iComponent);
        else
            error('Wall Condition not implemented')
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'fixedValue')
            cfdCorrectInletFixedValueBC(iPatch,theEquationName,iComponent);
        else
            error('Inlet Condition not implemented')
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalType,'outlet')
        if strcmp(theBCType,'fixedValue')
            cfdCorrectOutletFixedValueBC(iPatch,theEquationName,iComponent);
        elseif strcmp(theBCType,'fixedGradient')
            cfdCorrectOutletFixedGradientBC(iPatch,theEquationName,iComponent);
        elseif(strcmp(theBCType,'outlet'))
            cfdCorrectOutletOutletBC(iPatch,theEquationName,iComponent);
        elseif(strcmp(theBCType,'zeroGradient'))
            cfdCorrectOutletZeroGradientBC(iPatch,theEquationName,iComponent);
        else
            error([theBCType '<<<< Outlet Condition not implemented'])
        end
        %
        % SYMMETRY
        %
    elseif strcmp(thePhysicalType,'symmetry')
        cfdCorrectSymmetryBC(iPatch,theEquationName,iComponent);
        %
        % EMPTY
        %
    elseif strcmp(thePhysicalType,'empty')
        cfdCorrectEmptyBC(iPatch,theEquationName,iComponent);
        %
        % ERROR
        %
    elseif(strcmp(thePhysicalType,'patch'))
        if(strcmp(theBCType,'fixedValue'))
            cfdCorrectPatchSpecifiedValueBC(iPatch,theEquationName,iComponent);
        elseif(strcmp(theBCType,'fixedGradient'))
            
        elseif(strcmp(theBCType,'zeroGradient')||strcmp(theBCType,'Zero Gradient'))
            
        elseif(strcmp(theBCType,'outlet')|| strcmp(theBCType,'outlet'))
            
        else
            error([theBCType 'Wall bc not defined']);
        end
    else
        error('Symmetry Condition not implemented')
    end
    
end

end



%===================================================
% WALL-specifiedValue
%===================================================
function cfdCorrectWallFixedValueBC(iPatch,theEquationName,iComponent)
%

%
end

%===================================================
% WALL-zeroFlux
%===================================================
function cfdCorrectWallZeroGradientBC(iPatch,theEquationName,iComponent)
%
if(nargin==2)
    iComponent=1;
end
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi;
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

iOwners = [theMesh.faces(iBFaces).iOwner];
phi(iBElements,iComponent) = phi(iOwners,iComponent);
% PatchFACES
%
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% WALL-SPECIFIED FLUX
%===================================================
function cfdCorrectWallFixedGradientBC(iPatch,theEquationName,iComponent)
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi(:,iComponent);
%
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
%
% compute phiFlux
theEquation = cfdGetModel(theEquationName);
formula = theEquation.bcs{iPatch}.value;
theLocale = ['BPatch' num2str(iPatch)];
phiFlux(1:numberOfBFaces,iComponent)  = cfdComputeFormulaAtLocale(formula,theLocale);
% get gamma
theDiffusionTerm = cfdGetTermInEquation(theEquationName,'Diffusion');
theGammaField = cfdGetMeshField(theDiffusionTerm.coefficientName);
gamma = theGammaField.phi(iBElements,iComponent);
% get other info
gDiff = [theMesh.faces(iBFaces).gDiff]';
area = [theMesh.faces(iBFaces).area]';
%
iOwners = [theMesh.faces(iBFaces).iOwner];
%update boundary phi
phi(iBElements,iComponent) = (-phiFlux.*area - gamma.*gDiff.*phi(iOwners,iComponent))./(gamma.*gDiff);

%
% update meshfield
theMeshField.phi(:,iComponent) = phi;
cfdSetMeshField(theMeshField);

end

%===================================================
% WALL-HYBRID
%===================================================
function cfdCorrectWallHybridBC(iPatch,theEquationName)
%
theMesh = cfdGetMesh;
fvmFaces = theMesh.faces;
% boundary phi values
theGammaField = cfdGetModel(theScalarField.gamma);
gamma_b = theGammaField.phiPatches{iPatch};
%
phiFlux_b = theScalarField.phiFluxPatches{iPatch};
phi = theScalarField.phi;
phi_b = theScalarField.phiPatches{iPatch};
% PatchFACES
TheNumberOfPatchFaces = theMesh.patchSizes(iPatch);
patchFacesIndices = find([fvmFaces(:).patchIndex] == iPatch);

end

%===================================================
% WALL-noslip
%===================================================
function cfdCorrectWallNoslipBC(iPatch,theEquationName,iComponent)
%

%
end


%===================================================
% WALL-slip
%===================================================
function cfdCorrectWallSlipBC(iPatch,theEquationName,iComponent)
%
%
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

theModelEquation = cfdGetModel(theEquationName);

if(strcmp(theModelEquation.type,'Vector'))
    %
    %
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    %
    % Get Vector Field at Boundary
    %
    theVectorField = cfdGetMeshField(theEquationName);
    phi = theVectorField.phi;    
    %
    % update field
    %
    Sb = [theMesh.faces(iBFaces).Sf]';
    normSb = cfdMagnitude(Sb);
    n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];

    phi_normal = dot(phi(iOwners,:)',n')';
    
    phi(iBElements,iComponent) = phi(iOwners,iComponent) - phi_normal .* n(:,iComponent);
    
    theVectorField.phi = phi;
    cfdSetMeshField(theVectorField);       
else
    theScalarField = cfdGetMeshField(theEquationName);
    phi = theScalarField.phi;    
    %
    % PatchFACES
    %
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    phi(iBElements) = phi(iOwners);    
    %
    % update meshfield
    %
    theScalarField.phi = phi;
    cfdSetMeshField(theScalarField);
end
%
end


%===================================================
% INLET-specifiedValue
%===================================================
function cfdCorrectInletFixedValueBC(iPatch,theEquationName,iComponent)

end

%===================================================
% INLET-INLET
%===================================================
function cfdCorrectInletInletBC(iPatch,theEquationName,iComponent)

theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi(:,iComponent);

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
%
iOwners = [theMesh.faces(iBFaces).iOwner];
phi(iBElements) = phi(iOwners);


theMeshField.phi(:,iComponent) = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% OUTLET-OUTLET
%===================================================
function cfdCorrectOutletOutletBC(iPatch,theEquationName,iComponent)
%
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi(:,iComponent);
%
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
%
iOwners = [theMesh.faces(iBFaces).iOwner];
%
% update boundary phi
%
phi(iBElements) = phi(iOwners);
%
% update meshfield
%
theMeshField.phi(:,iComponent) = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% OUTLET-zeroFlux
%===================================================
function cfdCorrectOutletZeroGradientBC(iPatch,theEquationName,iComponent)
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi(:,iComponent);
%
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
%
%
iOwners = [theMesh.faces(iBFaces).iOwner];
%update boundary phi
phi(iBElements) = phi(iOwners);
%
% update meshfield
%
theMeshField.phi(:,iComponent) = phi;
cfdSetMeshField(theMeshField);

end

%===================================================
% OUTLET-specifiedValue
%===================================================
function cfdCorrectOutletFixedValueBC(iPatch,theEquationName,iComponent)
%

end

%===================================================
% OUTLET-SPECIFIED FLUX
%===================================================
function cfdCorrectOutletSpecifiedFlux(iPatch,theEquationName,iComponent)
%
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi(:,iComponent);
%
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
%
% compute phiFlux
theEquation = cfdGetModel(theEquationName);
formula = theEquation.bcs{iPatch}.value;
theLocale = ['BPatch' num2str(iPatch)];
phiFlux(1:numberOfBFaces,1)  = cfdComputeFormulaAtLocale(formula,theLocale);
% get gamma
theDiffusionTerm = cfdGetTermInEquation(theEquationName,'Diffusion');
theGammaField = cfdGetMeshField(theDiffusionTerm.coefficientName);
gamma = theGammaField.phi(iBElements);
% get other info
gDiff = [theMesh.faces(iBFaces).gDiff]';
area = [theMesh.faces(iBFaces).area]';
%
iOwners = [theMesh.faces(iBFaces).iOwner];
%update boundary phi
phi(iBElements) = (-phiFlux.*area - gamma.*gDiff.*phi(iOwners))./(-gamma.*gDiff);

%
% update meshfield
theMeshField.phi(:,iComponent) = phi;
cfdSetMeshField(theMeshField);

end

%===================================================
% SYMMETRY-zeroFlux
%===================================================
function cfdCorrectSymmetryBC(iPatch,theEquationName,iComponent)

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

theModelEquation = cfdGetModel(theEquationName);

if(strcmp(theModelEquation.type,'Vector'))
    %
    %
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    %
    % Get Velocity Field at Boundary
    %
    theVectorField = cfdGetMeshField(theEquationName);
    phi = theVectorField.phi;    
    %
    % update velocity
    %
    Sb = [theMesh.faces(iBFaces).Sf]';
    normSb = cfdMagnitude(Sb);
    n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];
    
    phi_normal = dot(phi(iOwners,:)',n')';
    
    phi(iBElements,iComponent) = phi(iOwners,iComponent) - phi_normal .* n(:,iComponent);
    
    theVectorField.phi = phi;
    cfdSetMeshField(theVectorField);    
    
else
    theScalarField = cfdGetMeshField(theEquationName);
    phi = theScalarField.phi;    
    %
    % PatchFACES
    %
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    phi(iBElements) = phi(iOwners);    
    %
    % update meshfield
    %
    theScalarField.phi = phi;
    cfdSetMeshField(theScalarField);
end
%
end


%===================================================
% EMPTY-zeroFlux
%===================================================
function cfdCorrectEmptyBC(iPatch,theEquationName,iComponent)

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

theModelEquation = cfdGetModel(theEquationName);

if(strcmp(theModelEquation.type,'Vector'))
    %
    %
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    %
    % Get Velocity Field at Boundary
    %
    theVectorField = cfdGetMeshField(theEquationName);
    phi = theVectorField.phi;    
    %
    % update velocity
    %
    Sb = [theMesh.faces(iBFaces).Sf]';
    normSb = cfdMagnitude(Sb);
    n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];

    phi_normal = dot(phi(iOwners,:)',n')';
    
    phi(iBElements,iComponent) = phi(iOwners,iComponent) - phi_normal .* n(:,iComponent);
    
    theVectorField.phi = phi;
    cfdSetMeshField(theVectorField);        
else
    theScalarField = cfdGetMeshField(theEquationName);
    phi = theScalarField.phi;    
    %
    % PatchFACES
    %
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    phi(iBElements) = phi(iOwners);    
    %
    % update meshfield
    %
    theScalarField.phi = phi;
    cfdSetMeshField(theScalarField);
end

end
%===================================================
% WALL-specifiedValue
%===================================================
function cfdCorrectPatchSpecifiedValueBC(iPatch,theEquationName,iComponent)
%

%
end

