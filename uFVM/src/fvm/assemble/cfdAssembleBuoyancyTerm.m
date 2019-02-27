function cfdAssembleBuoyancyTerm(theTerm, iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles the buoyancy term based on Boussinesq
%   approximation rho = rho0*(1 - beta*(T - TRef))
%--------------------------------------------------------------------------

% Get mesh info
theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;
iElements = 1:theNumberOfElements;
vol = [theMesh.elements.volume]';

% Get requierd fields and constants
theTemperatureField = cfdGetMeshField('T');
T = theTemperatureField.phi;

theTermFields = identifyFields(theTerm.formula);
if any(ismember(theTermFields, 'rho'))
    theDensityField = cfdGetMeshField('rho');
    rho = theDensityField.phi;
else
    rho = ones(theNumberOfElements, 1);
end

TRef = cfdGetConstant('TRef');
g = cfdGetConstant('g');
beta = cfdGetConstant('beta');

% Calculate fluxes
theFluxes.FLUXCE(iElements,1) =  0;
theFluxes.FLUXVE(iElements,1) =  0;
theFluxes.FLUXTE(iElements,1) = theTerm.sign * vol .* rho(iElements) .* g(iComponent) * beta .*(T(iElements) - TRef);

theFluxes.FLUXCEOLD(iElements,1) = 0;
theFluxes.FLUXTEOLD(iElements,1) = 0;

% Assemble
cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);