function [rmsResidual, maxResidual] = cfdAssembleEquation(theEquationName,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles equation
%--------------------------------------------------------------------------

if(nargin==1)
    iComponent=1;
end
theEquation = cfdGetModel(theEquationName);
theNumberOfTerms = length(theEquation.terms);
if(theNumberOfTerms>0)    
    % Initialize algebraic coefficients of the equation (ac, anb, bc, ...)
    theMesh = cfdGetMesh;
    theCoefficients = cfdSetupCoefficients(theMesh.cconn, theMesh.csize);
    cfdSetCoefficients(theCoefficients);
    
    % Assemble Equation Terms
    [rmsResidual, maxResidual] = cfdAssembleEquationTerms(theEquationName,iComponent);
    
    % Post Assemble Equation
    cfdPostAssembleEquation(theEquationName,iComponent);
end


