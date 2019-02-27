function theFluxes = cfdAssembleExplicitStressTerm(theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles the term in the momentum equation which
%   involves the second part of the shear stress term
%--------------------------------------------------------------------------

% Get fields
theVelocityField = cfdGetMeshField('U');

dUidx = theVelocityField.phiGradient(:,1,iComponent);
dUidx_f = cfdInterpolateFromElementsToFaces('Average', dUidx);

dUidy = theVelocityField.phiGradient(:,2,iComponent);
dUidy_f = cfdInterpolateFromElementsToFaces('Average', dUidy);

dUidz = theVelocityField.phiGradient(:,3,iComponent);
dUidz_f = cfdInterpolateFromElementsToFaces('Average', dUidz);

theTermFields = identifyFields(theTerm.formula);
if any(ismember(theTermFields, 'mu'))    
    theViscosityField = cfdGetMeshField('mu');
    mu = theViscosityField.phi;
    mu_f = cfdInterpolateFromElementsToFaces('Average', mu);
    gamma_f = mu_f;
elseif any(ismember(theTermFields, 'nu'))
    theViscosityField = cfdGetMeshField('nu');
    nu = theViscosityField.phi;
    nu_f = cfdInterpolateFromElementsToFaces('Average', nu);    
    gamma_f = nu_f;
end

theMesh = cfdGetMesh;
iFaces = 1:theMesh.numberOfInteriorFaces+theMesh.numberOfBFaces;
Sf = [theMesh.faces(iFaces).Sf]';

%---------------------------------------------------
% Assemble Over Interior Faces
%---------------------------------------------------
theFluxes.FLUXC1f(iFaces,1) =  0;
theFluxes.FLUXC2f(iFaces,1) =  0;
theFluxes.FLUXVf(iFaces,1) =  0;
theFluxes.FLUXTf(iFaces,1) =  theTerm.sign * gamma_f .* (dUidx_f .* Sf(:, 1) + dUidy_f .* Sf(:, 2) + dUidz_f .* Sf(:, 3));

cfdAssembleIntoGlobalMatrixFaceFluxes(theFluxes);

end


