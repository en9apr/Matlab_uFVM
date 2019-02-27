function cfdAssemblePressureGradientTerm(theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles the pressure gradient term
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;

% pressure gradient
thePressureField = cfdGetMeshField('p');
p_grad = thePressureField.phiGradient;

iElements = 1:theNumberOfElements;
%---------------------------------------------------
% Assemble Over Interior Faces
%---------------------------------------------------
%
if iComponent==1
    e = [1;0;0];
elseif iComponent==2
    e = [0;1;0];
elseif iComponent==3
    e = [0;0;1];
end

volume = [theMesh.elements.volume]';

theFluxes.FLUXCE(iElements,1) =  0;
theFluxes.FLUXVE(iElements,1) =  0;
theFluxes.FLUXTE(iElements,1) =  theTerm.sign * volume.*(p_grad(iElements,:)*e);

theFluxes.FLUXCEOLD(iElements,1) = 0;
theFluxes.FLUXTEOLD(iElements,1) = 0;

cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);
end


