function cfdInitializeMixtureMdotField
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
% get mesh info
theMesh = cfdGetMesh;

theMdotName = 'Mdot';
theMdotField = cfdGetMeshField(theMdotName,'Faces');
theVelocityField = cfdGetMeshField('Velocity');

theDensityField = cfdGetMeshField('Density');
% ////////////////////////////////////////
vel = theVelocityField.phi;
density = theDensityField.phi;
%
% interpolate fields to faces
%
vel_f = cfdInterpolateFromElementsToFaces('Average',vel);
density_f = cfdInterpolateFromElementsToFaces('Average',density);
%
% Initialize mdot_f at interior faces
%
iFaces = 1:theMesh.numberOfFaces;
Sf = [theMesh.faces(iFaces).Sf]';
mdot_f = density_f.*dot(Sf',vel_f')';

theMdotField.phi = mdot_f;
cfdSetMeshField(theMdotField);
