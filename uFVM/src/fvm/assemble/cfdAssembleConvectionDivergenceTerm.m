function cfdAssembleConvectionDivergenceTerm(theEquationName,theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function claculates the convection divergence related fluxes
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;
%---------------------- Start Term Info ---------------------
%
% specifiy the Term Field
%
theTermFieldName = theTerm.variableName;
theTermField = cfdGetMeshField(theTermFieldName);
phi = theTermField.phi(:,iComponent);

theCoefficients = cfdGetCoefficients;

theEquationBaseName = cfdGetBaseName(theEquationName);

% Get effective divergence
effDiv = cfdComputeEffectiveDivergenceWithFluidTag(theTerm);

if strcmp(theEquationBaseName,'VF')
    theNumberOfElements = theMesh.numberOfElements;
    for iElement= 1:theNumberOfElements
        FLUXCE = max(effDiv(iElement),0.0) - effDiv(iElement);
        FLUXVE = 0;
        FLUXTE = 0;

        theCoefficients.ac(iElement) = theCoefficients.ac(iElement) + FLUXCE;
        theCoefficients.bc(iElement) = theCoefficients.bc(iElement) - FLUXTE;
    end    
else
    theNumberOfElements = theMesh.numberOfElements;
    for iElement=1:theNumberOfElements
        FLUXCE =   max(effDiv(iElement),0.0) - effDiv(iElement);
        FLUXVE = - max(effDiv(iElement),0.0) * phi(iElement);
        FLUXTE = FLUXCE * phi(iElement) + FLUXVE;

        theCoefficients.ac(iElement) = theCoefficients.ac(iElement) + FLUXCE;
        theCoefficients.bc(iElement) = theCoefficients.bc(iElement) - FLUXTE;        
    end    
end

cfdSetCoefficients(theCoefficients);

end