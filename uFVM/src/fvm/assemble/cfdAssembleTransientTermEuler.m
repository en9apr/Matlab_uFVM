function theFluxes = cfdAssembleTransientTermEuler(theEquationName, theTerm, iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function calculates the transient fluxes based on Euler's method
%   and assembles them to the global matrix
%--------------------------------------------------------------------------

% Get mesh and geometry info
theMesh = cfdGetMesh;
iElements = 1:theMesh.numberOfElements;
vol = [theMesh.elements(iElements).volume]';

deltaT = cfdGetDt;

% Initialize fluxes
theFluxes.FLUXCE(iElements, 1)    = 0;
theFluxes.FLUXCEOLD(iElements, 1) = 0;
theFluxes.FLUXTE(iElements, 1)    = 0;

% Get phi and phi_old
theEquationMeshField = cfdGetMeshField(theEquationName);
phi = theEquationMeshField.phi(iElements,iComponent);

theEquationMeshFieldOld = cfdGetMeshField(theEquationName, 'Elements', 'Step1');
phi_old = theEquationMeshFieldOld.phi(iElements,iComponent);

% Get the rho field
theRhoField = cfdGetMeshField(['rho_',theEquationName,'eq']);
theRhoFieldOld = cfdGetMeshField(['rho_',theEquationName,'eq'],'Elements','Step1');

rho = theRhoField.phi(iElements);
rho_old = theRhoFieldOld.phi(iElements);

% Calculate element fluxes
theFluxes.FLUXCE    =   theTerm.sign * vol .* rho / deltaT;
theFluxes.FLUXCEOLD = - theTerm.sign * vol .* rho_old / deltaT;
theFluxes.FLUXTE    =   theFluxes.FLUXCE .* phi' + theFluxes.FLUXCEOLD .* phi_old';

% Assemble into global matrix
cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);


