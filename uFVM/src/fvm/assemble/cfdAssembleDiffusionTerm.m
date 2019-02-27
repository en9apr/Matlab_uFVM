function cfdAssembleDiffusionTerm(theEquationName,theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles Diffusion term
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
%
% Assemble Over Interior Faces
%
theFluxes = cfdAssembleDiffusionTermInterior(theEquationName,theTerm,iComponent);
%
% Assemble Over Boundary Patches
%
theNumberOfPatches = theMesh.numberOfBoundaries;
theEquation = cfdGetModel(theEquationName);
%
for iPatch=1:theNumberOfPatches
    % find the Physical Type
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    %
    theBCType = theEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if(strcmp(thePhysicalType,'wall'))
        if strcmp(theBCType,'fixedValue')
            theFluxes =  cfdAssembleDiffusionTermWallFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'fixedGradient')
            theFluxes =  cfdAssembleDiffusionTermWallFixedGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'zeroGradient')
            theFluxes =  cfdAssembleDiffusionTermWallZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'slip')
            theFluxes =  cfdAssembleDiffusionTermWallSlipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'noSlip')
            theFluxes =  cfdAssembleDiffusionTermWallNoslipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'Hybrid')
            theFluxes =  cfdAssembleDiffusionTermWallHybridBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType '<<<<< Not Implemented']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'fixedValue')
            theFluxes = cfdAssembleDiffusionTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'totalPressure')
            theFluxes = cfdAssembleStressTermInletTotalPressureBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType '<<<<< Not Implemented']);
        end
        %
        % OUTLET
        %
    elseif(strcmp(thePhysicalType,'outlet'))
        if strcmp(theBCType,'outlet') || strcmp(theBCType,'zeroGradient')
            theFluxes =  cfdAssembleDiffusionTermOutletBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType '<<<<< Not Implemented']);
        end
        %
        % SYMMETRY
        %        
    elseif strcmp(thePhysicalType,'symmetry')       
        theFluxes = cfdAssembleDiffusionTermSymmetry(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        %
        % EMPTY
        %        
    elseif strcmp(thePhysicalType,'empty')        
        theFluxes = cfdAssembleDiffusionTermEmpty(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        %
        % ERROR
        %         
    else
        error([thePhysicalType '<<<<<<< Not Implemented']);
    end
end

cfdAssembleIntoGlobalMatrixFaceFluxes(theFluxes);

end



%===================================================
% INTERIOR Faces
%===================================================
function theFluxes = cfdAssembleDiffusionTermInterior(theEquationName,theTerm,iComponent)
%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;

iFaces = 1:numberOfInteriorFaces;

%
%---------------------- Start Term Info ---------------------
%
% specifiy the Term Field
%
theTermFieldName = theTerm.variableName;
theTermField = cfdGetMeshField(theTermFieldName);
phi = theTermField.phi(:,iComponent);
phiGradient = theTermField.phiGradient(:,:,iComponent);

grad_f = cfdInterpolateGradientsFromElementsToInteriorFaces('Average:Corrected',phiGradient,phi);
grad_f = grad_f(iFaces,:);

% Get the term coefficient
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma__f = theGammaField.phi(iFaces);

gDiff_f = [theMesh.faces(iFaces).gDiff]';
Tf = [theMesh.faces(iFaces).T]';

iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';

theFluxes.FLUXC1f(iFaces,1) = - theTerm.sign * gamma__f .* gDiff_f;
theFluxes.FLUXC2f(iFaces,1) =   theTerm.sign * gamma__f .* gDiff_f;
theFluxes.FLUXVf(iFaces,1)  =   theTerm.sign * gamma__f .* dot(grad_f(:,:)',Tf(:,:)')';
theFluxes.FLUXTf(iFaces,1)  =   theFluxes.FLUXC1f(iFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iFaces) .* phi(iNeighbours) + theFluxes.FLUXVf(iFaces);


end

%===================================================
% INLET - fixedValue BC
%===================================================
function theFluxes = cfdAssembleDiffusionTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

% Get mesh info
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
gDiff = [theMesh.faces(iBFaces).gDiff]';
Tf = [theMesh.faces(iBFaces).T]';
iOwners = [theMesh.faces(iBFaces).iOwner]';

% Get field info
theTermFieldName = theTerm.variableName;
theTermField = cfdGetMeshField(theTermFieldName);
phi = theTermField.phi(:, iComponent);
grad = theTermField.phiGradient(:, :, iComponent);

% Get the gamma field
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);
   
% Add to face fluxes
theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * gamma .* gDiff;
theFluxes.FLUXC2f(iBFaces) =   theTerm.sign * gamma .* gDiff;
theFluxes.FLUXVf(iBFaces)  =   theTerm.sign * gamma .*dot(grad(iBElements,:)',Tf(:,:)')';
theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iBFaces) .* phi(iBElements) + theFluxes.FLUXVf(iBFaces);

end

%===================================================
% INLET - totalPressure BC
%===================================================
function theFluxes = cfdAssembleStressTermInletTotalPressureBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
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

end


%===================================================
% OUTLET - - BC
%===================================================
function theFluxes =  cfdAssembleDiffusionTermOutletBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;


theFluxes.FLUXC1f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXCVf(iBFaces) = 0;
theFluxes.FLUXTf(iBFaces)  = 0;

end


%===================================================
% WALL - fixedValue BC
%===================================================
function theFluxes = cfdAssembleDiffusionTermWallFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

% Get mesh info
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

% Get the term field
theTermFieldName = theTerm.variableName;
theTermField = cfdGetMeshField(theTermFieldName);
phi = theTermField.phi(:, iComponent);
grad = theTermField.phiGradient(:, :, iComponent);

% Get the gamma field
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);

% Get additional geometric info
gDiff = [theMesh.faces(iBFaces).gDiff]';
Tf = [theMesh.faces(iBFaces).T]';
iOwners = [theMesh.faces(iBFaces).iOwner]';
   
% Add to face fluxes
theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * gamma .* gDiff;
theFluxes.FLUXC2f(iBFaces) =   theTerm.sign * gamma .* gDiff;
theFluxes.FLUXVf(iBFaces)  =   theTerm.sign * gamma .* dot(grad(iBElements,:)',Tf(:,:)')';
theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iBFaces) .* phi(iBElements) + theFluxes.FLUXVf(iBFaces);

end


%===================================================
% WALL - fixedGradient BC
%===================================================
function theFluxes =  cfdAssembleDiffusionTermWallFixedGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
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
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;
%
%---------------------- Start Term Info ---------------------
%
% specifiy the Term Field
%
theTermFieldName = theTerm.variableName;
theTermField = cfdGetMeshField(theTermFieldName);
phi = theTermField.phi(:, iComponent);
grad = theTermField.phiGradient(:, :, iComponent);

%
% specifiy the Term Coefficient Field
%
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);

%---------------------- End Term Info ----------------------%
%
%
area = [theMesh.faces(iBFaces).area]';
   
%////////////////////////////////////////////////////////
theEquationField = cfdGetModel(theEquationName);
formula = theEquationField.bcs{iPatch}.value;
theLocale = ['BPatch' num2str(iPatch)];

% Calculate heat flux
normGrad_b = cfdComputeFormulaAtLocale(formula, theLocale);
phiFlux = - gamma .* normGrad_b;

%////////////////////////////////////////////////////////
%////////////////////////////////////////////////////////
%
theFluxes.FLUXC1f(iBFaces,1) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXTf(iBFaces)  = theTerm.sign * phiFlux .* area;

end

%===================================================
% WALL - zeroGradient BC
%===================================================
function theFluxes =  cfdAssembleDiffusionTermWallZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
theFluxes.FLUXC1f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXVf(iBFaces) = 0;
theFluxes.FLUXTf(iBFaces) = 0;

%
end


 
%===================================================
% WALL - HYBRID BC
%===================================================
function theFluxes =  cfdAssembleDiffusionTermWallHybridBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
theScalarField = cfdGetModel(theEquationName);
%
%---------------------- Start Term Info ---------------------
%
% specifiy the Term Field
%
if(isempty(theTerm.variableName))
    theTermField = theScalarField;
else
    theTermField = cfdGetModel(theTerm.variableName);
end
 
hc = cfdGetConstant('ConvectionCoefficient');
%Qradiation = cfdGetConstant('RadiationFlux');
Ta= cfdGetConstant('SolairTemperature');
 
%
% specifiy the Term Coefficient Field
%
theTermCoefficientField = fvmEvaluateFormula(theTerm.coefficient);
%---------------------- End Term Info ----------------------%
a = theFluxes.LHS;
b = theFluxes.RHS;
%
fvmFaces = theMesh.faces;
% boundary phi values
gamma_b = theTermCoefficientField.phiPatches{iPatch};
phi = theTermField.phi;
% loop over patchFaces
numberOfPatchFaces = theMesh.patchSizes(iPatch);
patchFacesIndices = find([fvmFaces(:).patchIndex] == iPatch);
%
 
   
   %    ////////////////////////////////////////////////////////
%    ////////////////////////////////////////////////////////
if(strcmp(theScalarField.type,'Energy')==true || strcmp(theScalarField.name,'H')==true)
    %
    cpField = cfdGetModel('SpecificHeat');
    LH = cfdGetConstant('LatentHeat');
    cp_b = cpField.phiPatches{iPatch};
 
 
    for k=1:numberOfPatchFaces
        iBFace = patchFacesIndices(k);
        theFace = fvmFaces(iBFace);
        iElement1 = theFace.element1;
 
        H_b = cp_b(k)*phi_b(k) + Vf_b(k)*LH ; % transform temperature bc into and enthalpy bc
 
        gDiff = theFace.gDiff;
        FLUXCb1 =   gamma_b(k)*gDiff;
        FLUXCb2 = - gamma_b(k)*gDiff;
        FLUXVb = 0;
        FLUXTb = FLUXCb1*phi(iElement1) + FLUXCb2*H_b + FLUXVb;
 
        a(iElement1,iElement1) = a(iElement1,iElement1) + cLHS*FLUXCb1;
        b(iElement1) = b(iElement1) - cRHS*FLUXTb;
    end
else
    
%This is the bit we edited

   for k=1:numberOfPatchFaces
        iBFace = patchFacesIndices(k);
        theFace = fvmFaces(iBFace);
        iElement1 = theFace.element1;
 
        gDiff = theFace.gDiff;
        FLUXCb1 =   gamma_b(k)*gDiff* gamma_b(k)*gDiff/(gamma_b(k)*gDiff+hc);
        FLUXCb2 = - gamma_b(k)*gDiff*hc/(gamma_b(k)*gDiff+hc);
        FLUXVb = 0;
        FLUXTb = FLUXCb1*phi(iElement1) + FLUXCb2*Ta + FLUXVb;
 
        a(iElement1,iElement1) = a(iElement1,iElement1) + cLHS*FLUXCb1;
        b(iElement1) = b(iElement1) - cRHS*FLUXTb;
    end
    
end
%
theFluxes.FLUXC1f(iBFaces) = FLUXCb1;
theFluxes.FLUXC2f(iBFaces) = FLUXCb2;
theFluxes.FLUXTf(iBFaces)  = FLUXTb;
 
end



%===================================================
% WALL - slip
%===================================================
function theFluxes =  cfdAssembleDiffusionTermWallSlipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

% Get mesh info
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXVf(iBFaces) = 0;
theFluxes.FLUXTf(iBFaces)  = 0;

end


%===================================================
% WALL - Noslip
%===================================================
function theFluxes =  cfdAssembleDiffusionTermWallNoslipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

% Get mesh info
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;

theEquationModel = cfdGetModel(theEquationName);
if strcmp(theEquationModel.type, 'Scalar')
    theFluxes.FLUXC1f(iBFaces) = 0;
    theFluxes.FLUXC2f(iBFaces) = 0;
    theFluxes.FLUXVf(iBFaces) = 0;
    theFluxes.FLUXTf(iBFaces)  = 0;
else

    % Get the vector field
    theVectorField = cfdGetMeshField(theEquationName);
    
    ex = [1;0;0];
    ey = [0;1;0];
    ez = [0;0;1];
    
    % Get the gamma field of the current equation
    theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
    gamma = theGammaField.phi(iBFaces);
    
    %---------------------- End Term Info ----------------------%
    %
    %
    area = [theMesh.faces(iBFaces).area]';
    n = [theMesh.faces(iBFaces).Sf]'./[area area area];
    dnorm = [theMesh.faces(iBFaces).walldist]';
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    
    TM = cfdGetUCoef(iPatch);
    
    %////////////////////////////////////////////////////////
    phi_c = theVectorField.phi(iOwners,:);
    phi_b =  theVectorField.phi(iBElements,:);
    
    nx = n*ex;
    ny = n*ey;
    nz = n*ez;
    
    if iComponent==1
        ac = - theTerm.sign * gamma .* area ./ dnorm .* TM .* (1- dot(nx',nx')');
        
        bc = - theTerm.sign * gamma .* area ./ dnorm .* TM.* (phi_b(:,1).*(1 - dot(nx',nx')') + ...
            (phi_c(:,2) - phi_b(:,2)).*ny.*nx - ...
            (phi_c(:,3) - phi_b(:,3)).*nz.*nx);
    elseif iComponent==2
        ac = - theTerm.sign * gamma .* area ./ dnorm .* TM .* (1- dot(ny',ny')');
        
        bc = - theTerm.sign * gamma .* area ./ dnorm .* TM.* (phi_b(:,2).*(1 - dot(ny',ny')') + ...
            (phi_c(:,1) - phi_b(:,1)).*nx.*ny - ...
            (phi_c(:,3) - phi_b(:,3)).*nz.*ny);
    else
        ac = - theTerm.sign * gamma .* area ./ dnorm .* TM .* (1- dot(nz',nz')');
        
        bc = - theTerm.sign * gamma .* area ./ dnorm .* TM.* (phi_b(:,3).*(1 - dot(nz',nz')') + ...
            (phi_c(:,1) - phi_b(:,1)).*nx.*nz - ...
            (phi_c(:,2) - phi_b(:,2)).*ny.*nz);
    end
    
    theFluxes.FLUXC1f(iBFaces) = ac;
    theFluxes.FLUXC2f(iBFaces) = 0;
    theFluxes.FLUXVf(iBFaces)  = bc;
    theFluxes.FLUXTf(iBFaces)  = theFluxes.FLUXC1f(iBFaces) .* phi_c(:,iComponent) + theFluxes.FLUXC2f(iBFaces) .* phi_b(:,iComponent) + theFluxes.FLUXVf(iBFaces);
    
end

%
end

%===================================================
% Symmetry - - BC
%===================================================
function theFluxes = cfdAssembleDiffusionTermSymmetry(iPatch,theFluxes,theEquationName,theTerm,iComponent)

% Get mesh info
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;

% Get the mesh field
theEquationModel = cfdGetModel(theEquationName);
if strcmp(theEquationModel.type, 'Scalar')
    theFluxes.FLUXC1f(iBFaces) = 0;
    theFluxes.FLUXC2f(iBFaces) = 0;
    theFluxes.FLUXVf(iBFaces)  = 0;
    theFluxes.FLUXTf(iBFaces)  = 0;
else
    theVectorField = cfdGetMeshField(theEquationName);
    phi = theVectorField.phi;
    
    if(iComponent==1)
        e = [1;0;0];
        l = repmat([0,1,1],numberOfBFaces,1);
    elseif(iComponent==2)
        e = [0;1;0];
        l = repmat([1,0,1],numberOfBFaces,1);
    elseif(iComponent==3)
        e = [0;0;1];
        l = repmat([1,1,0],numberOfBFaces,1);
    end
    
    % get the gamma field
    theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
    gamma_b = theGammaField.phi(iBFaces);
    
    % Get additional geometric info
    area = [theMesh.faces(iBFaces).area]';
    n = [theMesh.faces(iBFaces).Sf]'./[area area area];
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    CN = [theMesh.faces(iBFaces).CN]';
    d = dot(CN',n')';
    
    % Add to face fluxes
    theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * 2 * gamma_b .* area .* n(:,iComponent).^2 ./ d;
    theFluxes.FLUXC2f(iBFaces) =   0;
    theFluxes.FLUXVf(iBFaces)  =   theTerm.sign * 2 * gamma_b .* area .* dot((phi(iOwners,:).*l)',(n.*l)')' .* n(:,iComponent) ./ d;
    theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* phi(iOwners,iComponent) + theFluxes.FLUXC2f(iBFaces) .* phi(iBElements,iComponent) + theFluxes.FLUXVf(iBFaces);
    
end

end

%===================================================
% Empty - - BC
%===================================================
function theFluxes = cfdAssembleDiffusionTermEmpty(iPatch,theFluxes,theEquationName,theTerm,iComponent)

% Get mesh info
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;

iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;
%
iElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iElementEnd = iElementStart+numberOfBFaces-1;
iBElements = iElementStart:iElementEnd;

% Get the mesh field
theEquationModel = cfdGetModel(theEquationName);
if strcmp(theEquationModel.type, 'Scalar')
    theFluxes.FLUXC1f(iBFaces) = 0;
    theFluxes.FLUXC2f(iBFaces) = 0;
    theFluxes.FLUXVf(iBFaces)  = 0;
    theFluxes.FLUXTf(iBFaces)  = 0;
else
    theVectorField = cfdGetMeshField(theEquationName);
    phi = theVectorField.phi;
    
    if(iComponent==1)
        e = [1;0;0];
        l = repmat([0,1,1],numberOfBFaces,1);
    elseif(iComponent==2)
        e = [0;1;0];
        l = repmat([1,0,1],numberOfBFaces,1);
    elseif(iComponent==3)
        e = [0;0;1];
        l = repmat([1,1,0],numberOfBFaces,1);
    end
    
    % get the gamma field
    theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
    gamma_b = theGammaField.phi(iBFaces);
    
    % Get additional geometric info
    area = [theMesh.faces(iBFaces).area]';
    n = [theMesh.faces(iBFaces).Sf]'./[area area area];
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    CN = [theMesh.faces(iBFaces).CN]';
    d = dot(CN',n')';
    
    % Add to face fluxes
    theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * 2*gamma_b .* area .* n(:,iComponent).^2 ./ d;
    theFluxes.FLUXC2f(iBFaces) =   0;
    theFluxes.FLUXVf(iBFaces)  =   theTerm.sign * 2*gamma_b .* area .* dot((phi(iOwners,:).*l)',(n.*l)')' .* n(:,iComponent) ./ d;
    theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* phi(iOwners,iComponent) + theFluxes.FLUXC2f(iBFaces) .* phi(iBElements,iComponent) + theFluxes.FLUXVf(iBFaces);

end

end

