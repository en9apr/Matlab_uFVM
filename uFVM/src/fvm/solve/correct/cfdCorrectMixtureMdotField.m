function cfdCorrectMixtureMdotField
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function corrects the phasic and mixture mdot_f fields
%--------------------------------------------------------------------------

theNumberOfFluids = cfdGetNumberOfFluids;

theMesh = cfdGetMesh;
theNumberOfFaces = theMesh.numberOfFaces;

theMdotField = cfdGetMeshField('mdot_f','Faces');
theMdotField.phi(1:theNumberOfFaces, 1) = 0;

for iFluid=1:theNumberOfFluids
    cfdCorrectMdotInterior(iFluid);
    cfdCorrectMdotBoundaryPatches(iFluid);
    
    % Get volume fractions
    theVFField = cfdGetMeshField(['alpha', num2str(iFluid)]);
    vf_f = cfdInterpolateFromElementsToFaces('Average', theVFField.phi);
    
    % Get mdot_f
    thePhasicMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)], 'Faces');
    
    % Calculate mixture mdot_f
    theMdotField.phi = theMdotField.phi + vf_f .* thePhasicMdotField.phi;
end
cfdSetMeshField(theMdotField);


end


%===================================================
% Correct MdotField @ INTERIOR
%===================================================
function cfdCorrectMdotInterior(iFluid)

theFluxes = cfdGetFluxes;
FLUXC1f = theFluxes.phasic_FLUXC1f(:,iFluid);
FLUXC2f = theFluxes.phasic_FLUXC2f(:,iFluid);

theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_f = theMdotField.phi;

% Get mesh info
theMesh = cfdGetMesh;
iFaces = 1:theMesh.numberOfInteriorFaces;
iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';

% PP Field
thePPField = cfdGetMeshField('PP');
pp = thePPField.phi;

% Correct
mdot_f(iFaces) = mdot_f(iFaces) + FLUXC1f(iFaces).*pp(iOwners) + FLUXC2f(iFaces).*pp(iNeighbours);

% Store mdot_f
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end


%===================================================
% Correct Boundary Patches
%===================================================
function cfdCorrectMdotBoundaryPatches(iFluid)

theMesh = cfdGetMesh;
thePressureEquation = cfdGetModel('p');
%
theNumberOfPatches = theMesh.numberOfPatches;
for iPatch=1:theNumberOfPatches
   
    theBoundary = theMesh.boundaries(iPatch);
    thePhysicalType = theBoundary.type;
    %
    theBCType = thePressureEquation.bcs{iPatch}.type;
    %
    % WALL
    %
    if strcmp(thePhysicalType,'wall')
        if strcmp(theBCType,'slip')
            cfdCorrectMdotWallSlipBC(iPatch,iFluid);    
        elseif strcmp(theBCType,'noSlip') || strcmp(theBCType,'zeroGradient')
            cfdCorrectMdotWallNoslipBC(iPatch,iFluid);              
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % INLET
    %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'inlet') || strcmp(theBCType,'zeroGradient')
            cfdCorrectMdotInletInletBC(iPatch,iFluid);    
        elseif strcmp(theBCType,'fixedValue') 
            cfdCorrectMdotInletFixedValueBC(iPatch,iFluid);    
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % OUTLET
    %
    elseif strcmp(thePhysicalType,'outlet') 
        if strcmp(theBCType,'outlet')
            cfdCorrectMdotOutletOutletBC(iPatch,iFluid); 
        elseif(strcmp(theBCType,'fixedValue')) 
            cfdCorrectMdotOutletFixedValueBC(iPatch,iFluid); 
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % SYMMETRY
    %
    elseif strcmp(thePhysicalType,'symmetry')
        if strcmp(theBCType,'symmetry')
          cfdCorrectMdotSymmetryBC(iPatch,iFluid);  
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % ERROR
    %
    elseif strcmp(thePhysicalType,'empty')
           cfdCorrectMdotEmptyBC(iPatch,iFluid);         
    else
        error([thePhysicalType ' physical condition not implemented']);
    end        
    
end


end





%===================================================
% WALL-slip
%===================================================
function cfdCorrectMdotWallSlipBC(iPatch,iFluid)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

% update mdot
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_f = theMdotField.phi;
mdot_f(iBFaces) = 0.0*mdot_f(iBFaces);
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end


%===================================================
% WALL-noslip
%===================================================
function cfdCorrectMdotWallNoslipBC(iPatch,iFluid)
%
end

%===================================================
% INLET-INLET
%===================================================
function cfdCorrectMdotInletInletBC(iPatch,iFluid)

applicationClass = cfdGetApplicationClass;
if strcmp(applicationClass, 'compressible')
    theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
    mdot_f = theMdotField.phi;
    
    theMesh = cfdGetMesh;
    %
    theBoundary = theMesh.boundaries(iPatch);
    numberOfBFaces = theBoundary.numberOfBFaces;
    
    %
    iFaceStart = theBoundary.startFace;
    iFaceEnd = iFaceStart+numberOfBFaces-1;
    iBFaces = iFaceStart:iFaceEnd;
    
    iOwners = [theMesh.faces(iBFaces).iOwner]';
    
    %
    thePPField = cfdGetMeshField('PP');
    pp_b = thePPField.phi(iOwners);
    
    % Get denisty field
    theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
    rho = theDensityField.phi;
    rho_b = rho(iOwners);

    % Get drhodp field
    theDrhodpField = cfdGetMeshField(['C_rho', num2str(iFluid)]);
    C_rho = theDrhodpField.phi;
    C_rho_b = C_rho(iOwners);
    
    % Correct bt adding compressible contribution
    mdot_f(iBFaces) = mdot_f(iBFaces) + (mdot_f(iBFaces) ./ rho_b) .* C_rho_b .* pp_b;
    
    theMdotField.phi = mdot_f;
    cfdSetMeshField(theMdotField);
end

end


%===================================================
% INLET-specifiedValue
%===================================================
function cfdCorrectMdotInletFixedValueBC(iPatch,iFluid)

end


%===================================================
% OUTLET-Outlet
%===================================================
function cfdCorrectMdotOutletOutletBC(iPatch,iFluid)
%
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_f = theMdotField.phi;

theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
rho_b = theDensityField.phi;

theMesh = cfdGetMesh;
%
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
% Get Velocity Field at Boundary
%
theVelocityField = cfdGetMeshField('U');
U_b = theVelocityField.phi(iBElements,:);

% ////////////////////////////////////////////////////////////
% update mdot
Sb = [theMesh.faces(iBFaces).Sf]';
mdot_f(iBFaces) = mdot_f(iBFaces) + rho_b(iBElements).*dot(U_b',Sb')';
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end

%===================================================
% OUTLET-specifiedValue
%===================================================
function cfdCorrectMdotOutletFixedValueBC(iPatch,iFluid)
%
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_f = theMdotField.phi;

theDensityField = cfdGetMeshField(['rho', num2str(iFluid)]);
density_b = theDensityField.phi;

theMesh = cfdGetMesh;
%
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
% Get Velocity Field at Boundary
%
theVelocityField = cfdGetMeshField('U');
U_b = theVelocityField.phi(iOwners,:);
%
% update velocity
%
Sb = [theMesh.faces(iBFaces).Sf]';

% update mdot
mdot_f(iBFaces) = density_b(iBElements) .* dot(U_b',Sb')';
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);


end

function cfdCorrectMdotSymmetryBC(iPatch,iFluid)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

% Update mdot_f
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_f = theMdotField.phi;
mdot_f(iBFaces) = 0.0*mdot_f(iBFaces);
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end

function cfdCorrectMdotEmptyBC(iPatch,iFluid)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

% uUpdate mdot
theMdotField = cfdGetMeshField(['mdot_f', num2str(iFluid)],'Faces');
mdot_f = theMdotField.phi;
mdot_f(iBFaces) = 0.0*mdot_f(iBFaces);
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end
