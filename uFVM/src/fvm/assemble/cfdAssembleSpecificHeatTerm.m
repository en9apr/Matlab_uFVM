function theFluxes = cfdAssembleSpecificHeatTerm(theTerm)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles the term in the energy equation which involves
%   the specific heat
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;
iElements = 1:theNumberOfElements;
vol = [theMesh.elements(iElements).volume]';

% Specific Heat Field
theSpecificHeatField = cfdGetMeshField('Cp');
Cp = theSpecificHeatField.phi;
CpGrad = cfdComputeGradientGauss0(Cp,theMesh);
CpGrad = CpGrad(iElements,:);

% Temperature field
theTemperatureField = cfdGetMeshField('T');
T = theTemperatureField.phi(iElements);

% Density field
theDensityField = cfdGetMeshField('rho');
rho = theDensityField.phi;
rho = rho(iElements);

% Velocity field
theVelocityField = cfdGetMeshField('U');
U = theVelocityField.phi(iElements,:);

isTransient = cfdIsTransient;
if isTransient
    theSpecificHeatFieldOld = cfdGetMeshField('Cp','Elements','Step1');
    Cp_old = theSpecificHeatFieldOld.phi;
    
    dt = cfdGetDt;
    transientTerm = (Cp(iElements) - Cp_old(iElements)) / dt;
else
    transientTerm = 0;
end


%---------------------------------------------------
% Assemble Over Elements
%---------------------------------------------------

theFluxes.FLUXCE(iElements,1) = 0;
theFluxes.FLUXVE(iElements,1) = 0;
theFluxes.FLUXTE(iElements,1) = theTerm.sign * rho .* T .* (transientTerm + U(:, 1) .* CpGrad(:, 1) + U(:, 2) .* CpGrad(:, 2) + U(:, 3) .* CpGrad(:, 3)) .* vol;
theFluxes.FLUXCEOLD(iElements,1) = 0;
theFluxes.FLUXTEOLD(iElements,1) = 0;

cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);

end