function cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles algebraic equation coefficients from the
%   contribution of the element fluxes of the current term of the equation
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;

theCoefficients = cfdGetCoefficients;

% Call coefficients
ac = theCoefficients.ac;
ac_old = theCoefficients.ac_old;
bc = theCoefficients.bc;

% Assemble element fluxes
for iElement = 1:numberOfElements
    ac(iElement)     = ac(iElement)     + theFluxes.FLUXCE(iElement);
    ac_old(iElement) = ac_old(iElement) + theFluxes.FLUXCEOLD(iElement);
    bc(iElement)     = bc(iElement)     - theFluxes.FLUXTE(iElement);
end

% Store updated coefficients
theCoefficients.ac = ac;
theCoefficients.ac_old = ac_old;
theCoefficients.bc = bc;

cfdSetCoefficients(theCoefficients);

end