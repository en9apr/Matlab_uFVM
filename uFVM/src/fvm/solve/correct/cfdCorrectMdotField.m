function cfdCorrectMdotField
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function corrects the mdot_f field
%--------------------------------------------------------------------------

% Correct mdot field at interior faces
cfdCorrectMdotInterior;

% Correct mdot field at boundary faces
cfdCorrectMdotBoundaryPatches;


end


%===================================================
% Correct MdotField @ INTERIOR
%===================================================
function cfdCorrectMdotInterior

theFluxes = cfdGetFluxes;
FLUXC1f = theFluxes.FLUXC1f;
FLUXC2f = theFluxes.FLUXC2f;

theMdotField = cfdGetMeshField('mdot_f','Faces');
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
function cfdCorrectMdotBoundaryPatches

% get mesh info
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
            cfdCorrectMdotWallSlipBC(iPatch);    
        elseif strcmp(theBCType,'noSlip') || strcmp(theBCType,'zeroGradient')
            cfdCorrectMdotWallNoslipBC(iPatch);              
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % INLET
    %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'inlet') || strcmp(theBCType,'zeroGradient')
            cfdCorrectMdotInletInletBC(iPatch);    
        elseif strcmp(theBCType,'fixedValue') 
            cfdCorrectMdotInletFixedValueBC(iPatch);    
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % OUTLET
    %
    elseif strcmp(thePhysicalType,'outlet') 
        if strcmp(theBCType,'outlet')
            cfdCorrectMdotOutletOutletBC(iPatch); 
        elseif(strcmp(theBCType,'fixedValue')) 
            cfdCorrectMdotOutletFixedValueBC(iPatch); 
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % SYMMETRY
    %
    elseif strcmp(thePhysicalType,'symmetry')
        if strcmp(theBCType,'symmetry')
          cfdCorrectMdotSymmetryBC(iPatch);  
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % ERROR
    %
    elseif strcmp(thePhysicalType,'empty')
           cfdCorrectMdotEmptyBC(iPatch);         
    else
        error([thePhysicalType ' physical condition not implemented']);
    end        
    
end


end





%===================================================
% WALL-slip
%===================================================
function cfdCorrectMdotWallSlipBC(iPatch)

% Get mesh info
theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

% update mdot
theMdotField = cfdGetMeshField('mdot_f','Faces');
mdot_f = theMdotField.phi;
mdot_f(iBFaces) = 0.0*mdot_f(iBFaces);
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end


%===================================================
% WALL-noslip
%===================================================
function cfdCorrectMdotWallNoslipBC(iPatch)
%
end

%===================================================
% INLET-INLET
%===================================================
function cfdCorrectMdotInletInletBC(iPatch)

applicationClass = cfdGetApplicationClass;
if strcmp(applicationClass, 'compressible')
    if nargin==1
        theMdotField = cfdGetMeshField('mdot_f','Faces');
        mdot_f = theMdotField.phi;
    else
        theMdotField = cfdGetMeshField(['mdot_f', num2str],'Faces');
        mdot_f = theMdotField.phi;
    end
    
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
    theDensityField = cfdGetMeshField('rho');
    rho = theDensityField.phi;
    rho_b = rho(iOwners);

    % Get drhodp field
    theDrhodpField = cfdGetMeshField('C_rho');
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
function cfdCorrectMdotInletFixedValueBC(iPatch)

end


%===================================================
% OUTLET-Outlet
%===================================================
function cfdCorrectMdotOutletOutletBC(iPatch)

theMdotField = cfdGetMeshField('mdot_f','Faces');
mdot_f = theMdotField.phi;

theDensityField = cfdGetMeshField('rho_Ueq');
rho_b = theDensityField.phi;

% Get mesh info
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
function cfdCorrectMdotOutletFixedValueBC(iPatch)

theMdotField = cfdGetMeshField('mdot_f','Faces');
mdot_f = theMdotField.phi;

theDensityField = cfdGetMeshField('rho_Ueq');
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

function cfdCorrectMdotSymmetryBC(iPatch)

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

% Update mdot
theMdotField = cfdGetMeshField('mdot_f','Faces');
mdot_f = theMdotField.phi;
mdot_f(iBFaces) = 0.0*mdot_f(iBFaces);
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end

function cfdCorrectMdotEmptyBC(iPatch)

% Get mesh info
theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

% uUpdate mdot
theMdotField = cfdGetMeshField('mdot_f','Faces');
mdot_f = theMdotField.phi;
mdot_f(iBFaces) = 0.0*mdot_f(iBFaces);
theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);

end
