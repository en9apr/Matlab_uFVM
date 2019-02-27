function cfdUpdateTemperatureField(theTemperatureModel,totalTemperatureValue)
%===================================================

%  written by the CFD Group @ AUB, Spring 2016
%===================================================

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;

theTemperatureEquation = cfdGetModel(theTemperatureModel);
theFluidTag = theTemperatureEquation.tag;

theTemperatureField = cfdGetMeshField(['Temperature' theFluidTag]);
T_previous = theTemperatureField.phi(1:numberOfElements);
To = eval(totalTemperatureValue)*ones(size(T_previous));

theVelocityField = cfdGetMeshField(['Velocity' theFluidTag]);
V = theVelocityField.phi(1:numberOfElements,:);
Vmag = cfdMagnitude(V);

theSpecificHeatField = cfdGetMeshField(['SpecificHeat' theFluidTag]);
Cp = theSpecificHeatField.phi(1:numberOfElements);

T = To - Vmag.^2./(2*Cp);

URF = theTemperatureEquation.urf; 

%---------------------------------
% Correct Interior Field
%---------------------------------
theTemperatureField.phi(1:numberOfElements) = T_previous + URF*(T - T_previous);
cfdSetMeshField(theTemperatureField);

%
%---------------------------------
% Correct Boundary patches
%---------------------------------
theNumberOfBoundaries = theMesh.numberOfBoundaries;
for iPatch=1:theNumberOfBoundaries
    % get Physical and Equation Boundary Conditions
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    theBCType = theTemperatureEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if(strcmp(thePhysicalType,'Wall') || strcmp(thePhysicalType,'wall'))
        if(strcmp(theBCType,'specifiedValue'))
            cfdCorrectWallSpecifiedValueBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'specifiedFlux'))
            cfdCorrectWallSpecifiedFluxBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'zeroFlux'))
            cfdCorrectWallZeroFluxBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'Hybrid'))
            cfdCorrectWallHybridBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'noSlip')||strcmp(theBCType,'noSlip')||strcmp(theBCType,'noSlip')||strcmp(theBCType,'no-slip'))
            cfdCorrectWallNoslipBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'slip')||strcmp(theBCType,'slip'))
            cfdCorrectWallslipBC(iPatch,theTemperatureModel);
        else
            error('Wall Condition not implemented')
        end
        %
        % INLET
        %
    elseif(strcmp(thePhysicalType,'Inlet') || strcmp(thePhysicalType,'inlet'))
        if(strcmp(theBCType,'specifiedValue'))
            cfdCorrectInletSpecifiedValueBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'Inlet'))
            cfdCorrectInletInletBC(iPatch,theTemperatureModel);
        else
            error('Inlet Condition not implemented')
        end
        %
        % OUTLET
        %
    elseif(strcmp(thePhysicalType,'outlet') || strcmp(thePhysicalType,'outlet'))
        if(strcmp(theBCType,'specifiedValue'))
            cfdCorrectOutletSpecifiedValueBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'specifiedFlux'))
            cfdCorrectOutletSpecifiedFluxBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'outlet'))
            cfdCorrectOutletOutletBC(iPatch,theTemperatureModel);
        elseif(strcmp(theBCType,'zeroFlux'))
            cfdCorrectOutletZeroFluxBC(iPatch,theTemperatureModel);
        else
            error([theBCType '<<<< Outlet Condition not implemented'])
        end
        %
        % SYMMETRY
        %
    elseif(strcmp(thePhysicalType,'Symmetry') || strcmp(thePhysicalType,'symmetry'))
        cfdCorrectSymmetryBC(iPatch,theTemperatureModel);
        %
        % EMPTY
        %
    elseif(strcmp(thePhysicalType,'empty'))
        cfdCorrectEmptyZeroFluxBC(iPatch,theTemperatureModel);
        %
        % ERROR
        %
    else
        error('Symmetry Condition not implemented')
    end
    
end

end



%===================================================
% WALL-specifiedValue
%===================================================
function cfdCorrectWallSpecifiedValueBC(iPatch,theEquationName)
%

%
end

%===================================================
% WALL-zeroFlux
%===================================================
function cfdCorrectWallZeroFluxBC(iPatch,theEquationName)
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
phi(iBElements) = phi(iOwners);
% PatchFACES
%
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% WALL-SPECIFIED FLUX
%===================================================
function cfdCorrectWallSpecifiedFluxBC(iPatch,theEquationName)
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi;
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
phiFlux(1:numberOfBFaces)  = cfdComputeFormulaAtLocale(formula,theLocale);
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
phi(iBElements) = (-phiFlux.*area - gamma.*gDiff.*phi(iOwners))./(gamma.*gDiff);

%
% update meshfield
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end

%===================================================
% WALL-HYBRID
%===================================================
function cfdCorrectWallHybridBC(iPatch,theEquationName)
%
% Incomplete
%

end

%===================================================
% WALL-noslip
%===================================================
function cfdCorrectWallNoslipBC(iPatch,theEquationName)
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
phi(iBElements) = phi(iOwners);
% PatchFACES
%
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);
end


%===================================================
% WALL-slip
%===================================================
function cfdCorrectWallslipBC(iPatch,theEquationName)
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
phi(iBElements) = phi(iOwners);
% PatchFACES
%
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% INLET-specifiedValue
%===================================================
function cfdCorrectInletSpecifiedValueBC(iPatch,theEquationName)
%


end

%===================================================
% INLET-INLET
%===================================================
function cfdCorrectInletInletBC(iPatch,theEquationName)
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi;
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
phi(iBElements) = phi(iOwners);


theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% OUTLET-OUTLET
%===================================================
function cfdCorrectOutletOutletBC(iPatch,theEquationName)
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi;
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
%update boundary phi
phi(iBElements) = phi(iOwners);

%
% update meshfield
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% OUTLET-zeroFlux
%===================================================
function cfdCorrectOutletZeroFluxBC(iPatch,theEquationName)
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi(:);
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
%update boundary phi
phi(iBElements) = phi(iOwners);

%
% update meshfield
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end

%===================================================
% OUTLET-specifiedValue
%===================================================
function cfdCorrectOutletSpecifiedValueBC(iPatch,theEquationName)
%

end

%===================================================
% OUTLET-SPECIFIED FLUX
%===================================================
function cfdCorrectOutletSpecifiedFluxBC(iPatch,theEquationName)
%
%
% get phi
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi(:);
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
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end


%===================================================
% SYMMETRY-zeroFlux
%===================================================
function cfdCorrectSymmetryBC(iPatch,theEquationName)
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
% get phi
%
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi;
%
% PatchFACES
%
iOwners = [theMesh.faces(iBFaces).iOwner]';
phi(iBElements) = phi(iOwners);
%
% update meshfield
%
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);
%
end






%===================================================
% EMPTY-zeroFlux
%===================================================
function cfdCorrectEmptyZeroFluxBC(iPatch,theEquationName)
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
% get phi
%
theMeshField = cfdGetMeshField(theEquationName);
phi = theMeshField.phi;
%
% PatchFACES
%
iOwners = [theMesh.faces(iBFaces).iOwner]';
phi(iBElements) = phi(iOwners);
%
% update meshfield
%
theMeshField.phi = phi;
cfdSetMeshField(theMeshField);

end
