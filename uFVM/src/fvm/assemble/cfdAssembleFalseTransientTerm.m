function cfdAssembleFalseTransientTerm(theEquationName,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles false transient term
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
iElements = 1:theMesh.numberOfElements;

theEquation = cfdGetModel(theEquationName);
inertialTermExist = false;
for iTerm=1:length(theEquation.terms)
    termName = theEquation.terms{iTerm}.name;
    if strcmp(termName, 'Convection') || strcmp(termName, 'Transient')
        inertialTermExist = true;
        break;
    end
end

if ~inertialTermExist
    return;
end

% get rho field
theRhoField = cfdGetMeshField(['rho_',theEquationName,'eq']);
theRhoFieldOld = cfdGetMeshField(['rho_',theEquationName,'eq'],'Elements','Step1');

rho = theRhoField.phi(iElements);
rho_old = theRhoFieldOld.phi(iElements);

fdt = cfdGetFalseDt;

%---------------------- End Term Info ----------------------
theEquationMeshField = cfdGetMeshField(theEquationName);
theEquationMeshFieldOld = cfdGetMeshField(theEquationName, 'Elements', 'Step1');
phi = theEquationMeshField.phi(iElements,iComponent);
phi_old = theEquationMeshFieldOld.phi(iElements,iComponent);

volumes = [theMesh.elements(iElements).volume]';
theFluxes.FLUXCE(iElements,1)    =   volumes .* rho / fdt;
theFluxes.FLUXCEOLD(iElements,1) = - volumes .* rho_old / fdt;
theFluxes.FLUXTE(iElements,1)    =   theFluxes.FLUXCE(iElements,1) .* phi;
theFluxes.FLUXTEOLD(iElements,1) =   theFluxes.FLUXCEOLD(iElements,1) .* phi_old;

cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);
