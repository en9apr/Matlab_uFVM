function cfdAssembleGradientTerm(theEquationName,theTerm,iComponent)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
if(nargin==2)
    iComponent=1;
end

theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;

% pressure gradient
theField = cfdGetMeshField(theTerm.variableName);
if isempty(theField)
    phi = cfdComputeFormulaAtLocale(theTerm.variableName, 'Elements', 'Scalar', iComponent)';
    grad = cfdComputeFieldGradient(phi, theTerm.scheme);
else
    grad = theField.phiGradient;
end

iElements = 1:theNumberOfElements;
%---------------------------------------------------
% Assemble Over Interior Faces
%---------------------------------------------------
%
if(iComponent==1)
    e = [1;0;0];
elseif(iComponent==2)
    e = [0;1;0];
elseif(iComponent==3)
    e = [0;0;1];
end

volume = [theMesh.elements.volume]';

theFluxes.FLUXCE(iElements,1) =  0;
theFluxes.FLUXVE(iElements,1) =  0; % non orthogonal term
theFluxes.FLUXTE(iElements,1) =  theTerm.sign * volume.*(grad(iElements,:)*e);

theFluxes.FLUXCEOLD(iElements,1) = 0;
theFluxes.FLUXTEOLD(iElements,1) = 0;

cfdAssembleIntoGlobalMatrixElementFluxes(theEquationName,theFluxes,iComponent);
end


