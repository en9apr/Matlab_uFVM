function cfdAssembleStressTerm(theEquationName,theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles Stress term
%--------------------------------------------------------------------------

if(nargin==2)
    iComponent = 1;
end

theMesh = cfdGetMesh;
%---------------------------------------------------
% Assemble Over Interior Faces
%---------------------------------------------------
%
% Assemble Over Interior Faces
theFluxes  = cfdAssembleStressTermInterior(theEquationName,theTerm,iComponent);
%
% Assemble Over Boundary Patches
%
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
        if strcmp(theBCType,'slip')
            theFluxes =  cfdAssembleStressTermWallSlipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'noSlip')
            theFluxes =  cfdAssembleStressTermWallNoslipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif(strcmp(theBCType,'fixedValue'))
            theFluxes =  cfdAssembleStressTermWallFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType ' bc is not implemented']);
        end
        %
        % INLET
        %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'fixedValue')
            theFluxes =  cfdAssembleStressTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        elseif strcmp(theBCType,'totalPressure')
            theFluxes =  cfdAssembleStressTermInletTotalPressureBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType ' bc is not implemented']);
        end
        %
        % OUTLET
        %
    elseif strcmp(thePhysicalType,'outlet')
        if strcmp(theBCType,'zeroGradient')
            theFluxes =  cfdAssembleStressTermOutletZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        else
            error([theBCType ' bc is not implemented']);
        end
        %
        % SYMMETRY
        %
    elseif strcmp(thePhysicalType,'symmetry')
        theFluxes = cfdAssembleStressTermSymmetry(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        %
        % EMPTY
        %
    elseif strcmp(thePhysicalType,'empty')
        theFluxes = cfdAssembleStressTermEmptyBC(iPatch,theFluxes,theEquationName,theTerm,iComponent);
        %
        % ERROR
        %
    else
        error([thePhysicalType ' physical condition is not implemented']);
    end    
end

cfdAssembleIntoGlobalMatrixFaceFluxes(theFluxes);

end




%===================================================
% INTERIOR
%===================================================
function theFluxes = cfdAssembleStressTermInterior(theEquationName,theTerm,iComponent)

%
theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
iFaces = 1:numberOfInteriorFaces;
%
%---------------------- Start Term Info ---------------------
%
% specifiy the Term Field
%
theVelocityField = cfdGetMeshField('U');
phi = theVelocityField.phi(:,iComponent);
phiGradient = theVelocityField.phiGradient(:,:,iComponent);

grad_f = cfdInterpolateGradientsFromElementsToInteriorFaces('Average:Corrected',phiGradient,phi);
grad_f = grad_f(iFaces,:);
%
% specify the Term Coefficient Field
%
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iFaces);
%
%
% Assemble Interior Faces
%
geoDiff_f = [theMesh.faces(iFaces).gDiff]';
Tf = [theMesh.faces(iFaces).T]';

iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';

theFluxes.FLUXC1f(iFaces,1) = - theTerm.sign * gamma .* geoDiff_f;
theFluxes.FLUXC2f(iFaces,1) =   theTerm.sign * gamma .* geoDiff_f;
theFluxes.FLUXVf(iFaces,1)  =   theTerm.sign * gamma .* dot(grad_f(:,:)',Tf(:,:)')';  % non orthogonal term

theFluxes.FLUXVf(iFaces,1) = theFluxes.FLUXVf(iFaces,1);
theFluxes.FLUXTf(iFaces,1)  =  theFluxes.FLUXC1f(iFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iFaces) .* phi(iNeighbours) + theFluxes.FLUXVf(iFaces);

end

%===================================================
% WALLcfdAssembleStressTerm
%===================================================
function theFluxes =  cfdAssembleStressTermWallNoslipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

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
%
% specifiy the Term Field
%
theVelocityField = cfdGetMeshField('U');

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

%
% specifiy the Term Coefficient Field
%
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);

%---------------------- End Term Info ----------------------%
%
%
Sf = [theMesh.faces(iBFaces).Sf]';
area = [theMesh.faces(iBFaces).area]';
n = [theMesh.faces(iBFaces).Sf]'./[area area area];
dnorm = [theMesh.faces(iBFaces).walldist]';
iOwners = [theMesh.faces(iBFaces).iOwner]';

TM = cfdGetUCoef(iPatch);

%////////////////////////////////////////////////////////
velc = theVelocityField.phi(iOwners,:);
vel_b =  theVelocityField.phi(iBElements,:);

nx = n*ex;
ny = n*ey;
nz = n*ez;

if iComponent==1
    ac = - theTerm.sign * gamma .* area ./ dnorm .* TM .* (1- dot(nx',nx')');
    
    bc = - theTerm.sign * gamma .* area ./ dnorm .* TM.* (vel_b(:,1).*(1 - dot(nx',nx')') + ...
        (velc(:,2) - vel_b(:,2)).*ny.*nx - ...
        (velc(:,3) - vel_b(:,3)).*nz.*nx);
elseif iComponent==2
    ac = - theTerm.sign * gamma .* area ./ dnorm .* TM .* (1- dot(ny',ny')');
    
    bc = - theTerm.sign * gamma .* area ./ dnorm .* TM.* (vel_b(:,2).*(1 - dot(ny',ny')') + ...
        (velc(:,1) - vel_b(:,1)).*nx.*ny - ...
        (velc(:,3) - vel_b(:,3)).*nz.*ny);
else
    ac = - theTerm.sign * gamma .* area ./ dnorm .* TM .* (1- dot(nz',nz')');
    
    bc = - theTerm.sign * gamma .* area ./ dnorm .* TM.* (vel_b(:,3).*(1 - dot(nz',nz')') + ...
        (velc(:,1) - vel_b(:,1)).*nx.*nz - ...
        (velc(:,2) - vel_b(:,2)).*ny.*nz);
end

%////////////////////////////////////////////////////////
theFluxes.FLUXC1f(iBFaces) = ac;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXVf(iBFaces)  = bc;
theFluxes.FLUXTf(iBFaces)  = theFluxes.FLUXC1f(iBFaces) .* velc(:,iComponent) + theFluxes.FLUXC2f(iBFaces) .* vel_b(:,iComponent) + theFluxes.FLUXVf(iBFaces);

end


function theFluxes =  cfdAssembleStressTermWallSlipBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
%
theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXVf(iBFaces)  = 0; % Add non-orthogonal  term
theFluxes.FLUXTf(iBFaces)  = 0;

end

%===================================================
% INLET
%===================================================
%
% Specified Velocity
%
function theFluxes =  cfdAssembleStressTermInletFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

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
iOwners = [theMesh.faces(iBFaces).iOwner]';
%
% specifiy the Term Field
%
theVelocityField = cfdGetMeshField('U');
velc = theVelocityField.phi(:,iComponent);
%
% specifiy the Term Coefficient Field
%
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);

%---------------------- End Term Info ----------------------%
%
%
geodiff = [theMesh.faces(iBFaces).gDiff]';
%
theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * gamma.*geodiff;
theFluxes.FLUXC2f(iBFaces) =   theTerm.sign * gamma.*geodiff;
theFluxes.FLUXVf(iBFaces)  =   0; % Add non-orthogonal  term
theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* velc(iOwners) + theFluxes.FLUXC2f(iBFaces) .* velc(iBElements) + theFluxes.FLUXVf(iBFaces);
%
end
%
% Inlet (Specified Total Pressure)
%
function theFluxes =  cfdAssembleStressTermInletTotalPressureBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
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
% OUTLET
%===================================================
function theFluxes =  cfdAssembleStressTermOutletZeroGradientBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)
theMesh = cfdGetMesh;

theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;

%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

theFluxes.FLUXC1f(iBFaces) = 0;
theFluxes.FLUXC2f(iBFaces) = 0;
theFluxes.FLUXVf(iBFaces)  = 0;
theFluxes.FLUXTf(iBFaces)  = 0;

end

%===================================================
% SYMMETRY
%===================================================
function theFluxes = cfdAssembleStressTermSymmetry(iPatch,theFluxes,theEquationName,theTerm,iComponent)

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
%
% specifiy the Term Field
%
theTermField = cfdGetMeshField('U');
vel = theTermField.phi;

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
%
% specifiy the Term Coefficient Field
%
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);

%---------------------- End Term Info ----------------------%
%

Sf = [theMesh.faces(iBFaces).Sf]';
area = [theMesh.faces(iBFaces).area]';
n = [theMesh.faces(iBFaces).Sf]'./[area area area];
iOwners = [theMesh.faces(iBFaces).iOwner]';
CN = [theMesh.faces(iBFaces).CN]';
d = dot(CN',n')';
%////////////////////////////////////////////////////////

%////////////////////////////////////////////////////////
theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * 2*gamma .* area .* n(:,iComponent).^2 ./ d;
theFluxes.FLUXC2f(iBFaces) =   0;
theFluxes.FLUXVf(iBFaces)  =   theTerm.sign * 2*gamma .* area .* dot((vel(iOwners,:).*l)',(n.*l)')' .* n(:,iComponent) ./ d;
theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* vel(iOwners,iComponent) + theFluxes.FLUXC2f(iBFaces) .* vel(iBElements,iComponent) + theFluxes.FLUXVf(iBFaces);

end


%===================================================
% EMPTY
%===================================================
function theFluxes =  cfdAssembleStressTermEmptyBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

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
%
% specifiy the Term Field
%
theTermField = cfdGetMeshField('U');
velc = theTermField.phi;

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
%
% specifiy the Term Coefficient Field
%
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);

%---------------------- End Term Info ----------------------%
%

Sf = [theMesh.faces(iBFaces).Sf]';
area = [theMesh.faces(iBFaces).area]';
n = [theMesh.faces(iBFaces).Sf]'./[area area area];
iOwners = [theMesh.faces(iBFaces).iOwner]';
CN = [theMesh.faces(iBFaces).CN]';
d = dot(CN',n')';


theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * 2*gamma .* area .* n(:,iComponent).^2 ./ d;
theFluxes.FLUXC2f(iBFaces) =   0;
theFluxes.FLUXVf(iBFaces)  =   2*gamma .* area .* dot((velc(iOwners,:).*l)',(n.*l)')' .* n(:,iComponent) ./ d;
theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* velc(iOwners,iComponent) + theFluxes.FLUXC2f(iBFaces) .* velc(iBElements,iComponent) + theFluxes.FLUXVf(iBFaces);

end

%===================================================
% Wall- specifiedValue
%===================================================
function  theFluxes =  cfdAssembleStressTermWallFixedValueBC(iPatch,theFluxes,theEquationName,theTerm,iComponent)

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
%
% specifiy the Term Field
%
theTermField = cfdGetMeshField('U');
phi = theTermField.phi(:,iComponent);
%
% specifiy the Term Coefficient Field
%
theGammaField = cfdGetMeshField(['gamma_',theEquationName,'eq'], 'Faces');
gamma = theGammaField.phi(iBFaces);

%---------------------- End Term Info ----------------------%
%

geoDiff = [theMesh.faces(iBFaces).gDiff]';
Tf = [theMesh.faces(iBFaces).T]';
area = [theMesh.faces(iBFaces).area]';
n = [theMesh.faces(iBFaces).Sf]'./[area area area];
iOwners = [theMesh.faces(iBFaces).iOwner]';


theFluxes.FLUXC1f(iBFaces) = - theTerm.sign * gamma .* geoDiff ;
theFluxes.FLUXC2f(iBFaces) =   theTerm.sign *gamma .* geoDiff ;
theFluxes.FLUXVf(iBFaces)  =   0;
theFluxes.FLUXTf(iBFaces)  =   theFluxes.FLUXC1f(iBFaces) .* phi(iOwners) + theFluxes.FLUXC2f(iBFaces) .* phi(iBElements) + theFluxes.FLUXVf(iBFaces);
end
