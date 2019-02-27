function cfdCreateMdotField
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function initializes mdot_f term
%--------------------------------------------------------------------------

% get mesh and field info
theMesh = cfdGetMesh;
theVelocityField = cfdGetMeshField('U');
theDensityField = cfdGetMeshField('rho');
if isempty(theDensityField)
    rho = ones(size(theVelocityField.phi(:,1)));
else
    rho = theDensityField.phi;
end
theMdotField = cfdGetMeshField('mdot_f', 'Faces');

if cfdIsTransient
    % Get required fields
    dt = cfdGetDt;
    U = theVelocityField.phi;
    
    % interpolate fields to faces
    U_f = cfdInterpolateFromElementsToFaces('Average', U);
    rho_f = cfdInterpolateFromElementsToFaces('Average', rho);
    
    % Initialize mdot_f at interior faces
    iFaces = 1:theMesh.numberOfFaces;
    Sf = [theMesh.faces(iFaces).Sf]';
    DeltaVol = [theMesh.faces(iFaces).DeltaVol]';
    mdot_f1 = rho_f .* dot(Sf',U_f')';
    mdot_f2 = rho_f .* DeltaVol / dt;
    mdot_f = mdot_f1 - mdot_f2;
    
    theMdotField.phi = mdot_f;
    cfdSetMeshField(theMdotField);
else
    U = theVelocityField.phi;
    
    % interpolate fields to faces
    U_f = cfdInterpolateFromElementsToFaces('Average',U);
    rho_f = cfdInterpolateFromElementsToFaces('Average',rho);
    
    % Initialize mdot_f at interior faces
    iFaces = 1:theMesh.numberOfFaces;
    Sf = [theMesh.faces(iFaces).Sf]';
    mdot_f = rho_f.*dot(Sf',U_f')';
    
    theMdotField.phi = mdot_f;
    cfdSetMeshField(theMdotField);
end