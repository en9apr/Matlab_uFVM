function theFluxes = cfdAssembleConvectionTerm(theEquationName,theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles convection term
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
%---------------------------------------------------
% Assemble Over Interior Faces
%---------------------------------------------------
theFluxes = cfdAssembleConvectionTermInterior(theEquationName,theTerm,iComponent);
%---------------------------------------------------
% Assemble  Over Boundary Patches
%---------------------------------------------------
theNumberOfPatches = theMesh.numberOfBoundaries;
theEquation = cfdGetModel(theEquationName);

for iPatch=1:theNumberOfPatches
    % find the Physical Type
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    %
    theBCType = theEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if strcmp(thePhysicalType,'wall')
        if strcmp(theBCType,'slip') || strcmp(theBCType,'noSlip')
            theFluxes =  cfdAssembleConvectionTermWallBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'fixedGradient') || strcmp(theBCType,'fixedValue')
            theFluxes =  cfdAssembleConvectionTermWallBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType '<<<< Not implemented']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'fixedValue')
            theFluxes =  cfdAssembleConvectionTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType '<<<< Not implemented']);
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalType,'outlet')
        if strcmp(theBCType,'zeroGradient') || strcmp(theBCType,'outlet')
            theFluxes =  cfdAssembleConvectionTermOutletZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType '<<<< Not implemented']);
        end
        %
        % Symmetry
        %
    elseif strcmp(thePhysicalType,'symmetry')
        theFluxes =  cfdAssembleConvectionTermSymmetryBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        %
        % Empty
        %
    elseif strcmp(thePhysicalType,'empty')
        theFluxes =  cfdAssembleConvectionTermEmptyBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        %
        % ERROR
        %
    else
        error([thePhysicalType '<<<< Not implemented']);
    end
end

% Correct for Deferred Correction
theScheme = theTerm.scheme;
if strcmp(theScheme,'SOU')
    theFluxes = cfdAssembleConvectionTermDCSOU(theEquationName,theTerm,theFluxes,iComponent);
elseif strcmp(theScheme,'CENTRAL')
    theFluxes = cfdAssembleConvectionTermDCCENTRAL(theEquationName,theTerm,theFluxes,iComponent);
elseif strcmp(theScheme,'QUICK')
    theFluxes = cfdAssembleConvectionTermDCQUICK(theEquationName,theTerm,theFluxes,iComponent);
elseif strcmp(theScheme,'SMART')
    theFluxes = cfdAssembleConvectionTermDCSMART(theEquationName,theTerm,theFluxes,iComponent);
elseif strcmp(theScheme,'STOIC')
    theFluxes = cfdAssembleConvectionTermDCSTOIC(theEquationName,theTerm,theFluxes,iComponent);
elseif strcmp(theScheme,'MINMOD')
    theFluxes = cfdAssembleConvectionTermDCMINMOD(theEquationName,theTerm,theFluxes,iComponent);
end

% Assemble into global matrix
cfdAssembleIntoGlobalMatrixFaceFluxes(theFluxes);

% Assemble convection divergence term
cfdAssembleConvectionDivergenceTerm(theEquationName,theTerm,iComponent);

end


%===================================================
% Assemble Interior Faces
%===================================================
function theFluxes = cfdAssembleConvectionTermInterior(theEquationName,theTerm,iComponent)

% Get mesh info
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

iFaces = 1:numberOfInteriorFaces;
iOwners = [theMesh.faces(iFaces).iOwner];
iNeighbours = [theMesh.faces(iFaces).iNeighbour];

%
%---------------------- Start Term Info ---------------------
%
% Get the term variable
%
theVariableName = theTerm.variableName;
theVariableField = cfdGetMeshField(theVariableName);
phi = theVariableField.phi(:,iComponent);

% Get the psi field term coefficient. In case of velocity equation, the
% term coefficient is simply the mdot term
if strcmp(theEquationName, 'U')
    thePsiField = cfdGetMeshField('mdot_f', 'Faces');
else
    thePsiField = cfdGetMeshField(['psi_',theEquationName,'eq'], 'Faces');
end

psi = thePsiField.phi(iFaces);

theFluxes.FLUXC1f(iFaces,1) =   theTerm.sign * max(psi,0);
theFluxes.FLUXC2f(iFaces,1) = - theTerm.sign * max(-psi,0);
theFluxes.FLUXVf(iFaces,1) =    0;
theFluxes.FLUXTf(iFaces,1) = theFluxes.FLUXC1f(iFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iFaces) .* phi(iNeighbours) + theFluxes.FLUXVf(iFaces);

end



%===================================================
% INLET BC
%===================================================
function theFluxes = cfdAssembleConvectionTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
%---------------------- Start Term Info ---------------------
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
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
iOwners = [theMesh.faces(iBFaces).iOwner];
%
%---------------------- Start Term Info ---------------------
%
% specifiy the Term Field
%
theVariableName = theTerm.variableName;
theVariableField = cfdGetMeshField(theVariableName);
phi = theVariableField.phi(:,iComponent);

% Get the term coefficient
if strcmp(theEquationName, 'U')
    thePsiField = cfdGetMeshField('mdot_f', 'Faces');
else
    thePsiField = cfdGetMeshField(['psi_',theEquationName,'eq'], 'Faces');
end

psi = thePsiField.phi(iBFaces);

theFluxes.FLUXC1f(iBFaces,1) =   theTerm.sign * max(psi,0);
theFluxes.FLUXC2f(iBFaces,1) = - theTerm.sign * max(-psi,0);
theFluxes.FLUXVf(iBFaces,1) = 0;
theFluxes.FLUXTf(iBFaces,1)  = theFluxes.FLUXC1f(iBFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iBFaces) .* phi(iBElements) + theFluxes.FLUXVf(iBFaces);

end

%===================================================
% OUTLET BC
%===================================================
function theFluxes =  cfdAssembleConvectionTermOutletZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
%---------------------- Start Term Info ---------------------
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
iOwners = [theMesh.faces(iBFaces).iOwner];
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
%---------------------- Start Term Info ---------------------
%
% specifiy the Term Field
%
theVariableName = theTerm.variableName;
theVariableField = cfdGetMeshField(theVariableName);
phi = theVariableField.phi(:,iComponent);

% Get the term coefficient
if strcmp(theEquationName, 'U')
    thePsiField = cfdGetMeshField('mdot_f', 'Faces');
else
    thePsiField = cfdGetMeshField(['psi_',theEquationName,'eq'], 'Faces');
end

psi = thePsiField.phi(iBFaces);

theFluxes.FLUXC1f(iBFaces,1) =   theTerm.sign * max(psi,0);
theFluxes.FLUXC2f(iBFaces,1) = - theTerm.sign * max(-psi,0);
theFluxes.FLUXVf(iBFaces,1) = 0;
theFluxes.FLUXTf(iBFaces,1)  = theFluxes.FLUXC1f(iBFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iBFaces) .* phi(iBElements) + theFluxes.FLUXVf(iBFaces);

end


%===================================================
% WALL BC
%===================================================
function theFluxes =  cfdAssembleConvectionTermWallBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
%
%---------------------- Start Term Info ---------------------
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces,1) = 0;
theFluxes.FLUXC2f(iBFaces,1) = 0;
theFluxes.FLUXVf(iBFaces,1) = 0;
theFluxes.FLUXTf(iBFaces,1) = 0;

end


%===================================================
% EMPTY BC
%===================================================
function theFluxes =  cfdAssembleConvectionTermEmptyBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
%
%---------------------- Start Term Info ---------------------
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXVf(iBFaces) = 0;
theFluxes.FLUXTf(iBFaces) = 0;

end



%===================================================
% SYMMETRY BC
%===================================================
function theFluxes =  cfdAssembleConvectionTermSymmetryBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
%
%---------------------- Start Term Info ---------------------
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXVf(iBFaces)  = 0;
theFluxes.FLUXTf(iBFaces)  = 0;

end




