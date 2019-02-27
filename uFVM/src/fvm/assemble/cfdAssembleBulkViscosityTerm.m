function theFluxes = cfdAssembleBulkViscosityTerm(theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles the term in the compressible momentum equation
%   involving the bulk viscosity
%--------------------------------------------------------------------------

% Get fields
theVelocityField = cfdGetMeshField('U');

dudx = theVelocityField.phiGradient(:,1,1);
dudx_f = cfdInterpolateFromElementsToFaces('Average', dudx);

dvdy = theVelocityField.phiGradient(:,2,2);
dvdy_f = cfdInterpolateFromElementsToFaces('Average', dvdy);

dwdz = theVelocityField.phiGradient(:,3,3);
dwdz_f = cfdInterpolateFromElementsToFaces('Average', dwdz);

theViscosityField = cfdGetMeshField('mu');
mu = theViscosityField.phi;
mu_f = cfdInterpolateFromElementsToFaces('Average', mu);

theMesh = cfdGetMesh;
iFaces = 1:theMesh.numberOfInteriorFaces+theMesh.numberOfBFaces;
Sf = [theMesh.faces(iFaces).Sf]';

%---------------------------------------------------
% Assemble Over Interior Faces
%---------------------------------------------------
theFluxes.FLUXC1f(iFaces,1) =  0;
theFluxes.FLUXC2f(iFaces,1) =  0;
theFluxes.FLUXVf(iFaces,1) =  0;
theFluxes.FLUXTf(iFaces,1) =  theTerm.sign * 2/3 * mu_f .* (dudx_f + dvdy_f + dwdz_f) .* Sf(:, iComponent);

cfdAssembleIntoGlobalMatrixFaceFluxes(theFluxes);

end


