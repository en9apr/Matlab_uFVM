function theFluxes = cfdAssembleSubstantialDerivativeTerm(theTerm)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles the term in the energy equation which involves
%   the substantial derivarive of the pressure
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;
iElements = 1:theNumberOfElements;
vol = [theMesh.elements(iElements).volume]';

% Pressure field
thePressureField = cfdGetMeshField('p');
p = thePressureField.phi(iElements);
pGrad = thePressureField.phiGradient(iElements, :);

% Velocity field
theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi(iElements,:);

isTransient = cfdIsTransient;
if isTransient
    thePressureFieldOld = cfdGetMeshField('p','Elements','Step1');
    p_old = thePressureFieldOld.phi;
    
    dt = cfdGetDt;
    transientTerm = (p(iElements) - p_old(iElements)) / dt;
else
    transientTerm = 0;
end


%---------------------------------------------------
% Assemble Over Elements
%---------------------------------------------------

theFluxes.FLUXCE(iElements,1) = 0;
theFluxes.FLUXVE(iElements,1) = 0;
theFluxes.FLUXTE(iElements,1) = theTerm.sign * (transientTerm + U(:, 1) .* pGrad(:, 1) + U(:, 2) .* pGrad(:, 2) + U(:, 3) .* pGrad(:, 3)) .* vol;
theFluxes.FLUXCEOLD(iElements,1) = 0;
theFluxes.FLUXTEOLD(iElements,1) = 0;

cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);


end
