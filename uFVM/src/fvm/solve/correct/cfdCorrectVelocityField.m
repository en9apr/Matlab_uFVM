function cfdCorrectVelocityField
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function corrects the equations
%--------------------------------------------------------------------------

%------------------------------------
% Correct VelocityField @ Interior
%------------------------------------
cfdCorrectVelocityInterior;

%===================================================
% Correct Boundary Patches
%===================================================
cfdCorrectVelocityBoundaryPatches;

end

%------------------------------------
% Correct VelocityField @ Interior
%------------------------------------
function cfdCorrectVelocityInterior

theMesh = cfdGetMesh;

iElements = 1:theMesh.numberOfElements;

theD1Field = cfdGetMeshField('DU1');
theD2Field = cfdGetMeshField('DU2');
theD3Field = cfdGetMeshField('DU3');
DU1 = theD1Field.phi(iElements);
DU2 = theD2Field.phi(iElements);
DU3 = theD3Field.phi(iElements);
%  
% PP Field
%
thePPField = cfdGetMeshField('PP');
ppGrad = thePPField.phiGradient;
%
DUPPGRAD = [DU1.*ppGrad(iElements,1),DU2.*ppGrad(iElements,2),DU3.*ppGrad(iElements,3)];
%
theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi;
%
U(iElements,:) = U(iElements,:) - DUPPGRAD(iElements,:);

%
theVelocityField.phi = U;
cfdSetMeshField(theVelocityField);

end

function cfdCorrectVelocityBoundaryPatches
%
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
            cfdCorrectMdotAndVelocityWallSlipBC(iPatch);    
        elseif strcmp(theBCType,'noSlip') || strcmp(theBCType,'zeroGradient')
            cfdCorrectMdotAndVelocityWallNoslipBC(iPatch);              
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % INLET
    %
    elseif strcmp(thePhysicalType,'inlet')
        if strcmp(theBCType,'inlet') || strcmp(theBCType,'zeroGradient')
            cfdCorrectMdotAndVelocityInletInletBC(iPatch);    
        elseif strcmp(theBCType,'fixedValue') 
            cfdCorrectMdotAndVelocityInletFixedValueBC(iPatch);    
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % OUTLET
    %
    elseif strcmp(thePhysicalType,'outlet') 
        if strcmp(theBCType,'outlet')
            cfdCorrectMdotAndVelocityOutletOutletBC(iPatch); 
        elseif(strcmp(theBCType,'fixedValue')) 
            cfdCorrectMdotAndVelocityOutletFixedValueBC(iPatch); 
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % SYMMETRY
    %
    elseif strcmp(thePhysicalType,'symmetry')
        if strcmp(theBCType,'symmetry')
          cfdCorrectMdotAndVelocitySymmetryBC(iPatch);  
        else
            error([theBCType ' boundary condition not implemented']);
        end
    %
    % ERROR
    %
    elseif strcmp(thePhysicalType,'empty')
           cfdCorrectMdotAndVelocityEmptyBC(iPatch);  
       
    else
        error([thePhysicalType ' physical condition not implemented']);
    end        
    
end


end





%===================================================
% WALL-slip
%===================================================
function cfdCorrectMdotAndVelocityWallSlipBC(iPatch)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

iBElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iBElementEnd = iBElementStart+numberOfBFaces-1;
iBElements = iBElementStart:iBElementEnd;

iOwners = [theMesh.faces(iBFaces).iOwner];

theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi;

% Update velocity
Sb = [theMesh.faces(iBFaces).Sf]';
normSb = cfdMagnitude(Sb);
n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];

U_normal_mag = dot(U(iOwners,:)',n')';
U_normal = [U_normal_mag .* n(:,1),U_normal_mag .* n(:,2),U_normal_mag .* n(:,3)];

U(iBElements,:) = U(iOwners,:) - U_normal;

theVectorField.phi = U;
cfdSetMeshField(theVectorField);

end


%===================================================
% WALL-noslip
%===================================================
function cfdCorrectMdotAndVelocityWallNoslipBC(iPatch)
%
end

%===================================================
% INLET-INLET
%===================================================
function cfdCorrectMdotAndVelocityInletInletBC(iPatch)

end


%===================================================
% INLET-specifiedValue
%===================================================
function cfdCorrectMdotAndVelocityInletFixedValueBC(iPatch)
%
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
    
% get Fluxes 
theFluxes = cfdGetFluxes;
FLUXC1f = theFluxes.FLUXC1f;
FLUXC2f = theFluxes.FLUXC2f;

thePPField = cfdGetMeshField('PP');
pp = thePPField.phi;
%
% Get Velocity Field at Boundary
%
theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi;

theDensityField = cfdGetMeshField('rho_Ueq');
density_b = theDensityField.phi(iBElements);

% Update velcity
Sb = [theMesh.faces(iBFaces).Sf]';
normSb = norm(Sb);
ev = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];
dv_mag = (FLUXC1f(iBFaces).*pp(iOwners)./density_b)./dot(ev',Sb')'; 
%
U(iBElements,:) = U(iBElements,:) + [dv_mag.*ev(:,1), dv_mag.*ev(:,2), dv_mag.*ev(:,3)];
theVelocityField.phi = U;
cfdSetMeshField(theVelocityField);
%
end


%===================================================
% OUTLET-Outlet
%===================================================
function cfdCorrectMdotAndVelocityOutletOutletBC(iPatch)

end

%===================================================
% OUTLET-specifiedValue
%===================================================
function cfdCorrectMdotAndVelocityOutletFixedValueBC(iPatch)
%
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
iOwners = [theMesh.faces(iBFaces).iOwner]';    
%
% Get Velocity Field at Boundary
%
theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi;
U_b = U(iOwners,:);

theVelocityField.phi(iBFaces,:) = U_b;
cfdSetMeshField(theVelocityField);

end

function cfdCorrectMdotAndVelocitySymmetryBC(iPatch)
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

iBElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iBElementEnd = iBElementStart+numberOfBFaces-1;
iBElements = iBElementStart:iBElementEnd;

iOwners = [theMesh.faces(iBFaces).iOwner];

theVectorField = cfdGetMeshField('U');
U = theVectorField.phi;

% Update velocity
Sb = [theMesh.faces(iBFaces).Sf]';
normSb = cfdMagnitude(Sb);
n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];

U_normal_mag = dot(U(iOwners,:)',n')';
U_normal = [U_normal_mag .* n(:,1),U_normal_mag .* n(:,2),U_normal_mag .* n(:,3)];

U(iBElements,:) = U(iOwners,:) - U_normal;

theVectorField.phi = U;
cfdSetMeshField(theVectorField);

end

function cfdCorrectMdotAndVelocityEmptyBC(iPatch)

theMesh = cfdGetMesh;
theBoundary = theMesh.boundaries(iPatch);
numberOfElements = theMesh.numberOfElements;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfBFaces = theBoundary.numberOfBFaces;
%
iFaceStart = theBoundary.startFace;
iFaceEnd = iFaceStart+numberOfBFaces-1;
iBFaces = iFaceStart:iFaceEnd;

iBElementStart = numberOfElements+iFaceStart-numberOfInteriorFaces;
iBElementEnd = iBElementStart+numberOfBFaces-1;
iBElements = iBElementStart:iBElementEnd;

iOwners = [theMesh.faces(iBFaces).iOwner];

theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi;

% Update velocity
Sb = [theMesh.faces(iBFaces).Sf]';
normSb = cfdMagnitude(Sb);
n = [Sb(:,1)./normSb, Sb(:,2)./normSb, Sb(:,3)./normSb];

U_normal_mag = dot(U(iOwners,:)',n')';
U_normal = [U_normal_mag .* n(:,1),U_normal_mag .* n(:,2),U_normal_mag .* n(:,3)];

U(iBElements,:) = U(iOwners,:) - U_normal;

theVelocityField.phi = U;
cfdSetMeshField(theVelocityField);
end
