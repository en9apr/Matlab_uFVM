function cfdTransientUpdate
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function stores current time step fields to previous one
%--------------------------------------------------------------------------

% Update Scalar Fields
%
theEquationNames = cfdGetEquationNames;
theNumberOfEquations = length(theEquationNames);
for iEquation = 1:theNumberOfEquations
    theEquationName = theEquationNames{iEquation};
    cfdUpdateTransientEquation(theEquationName);
    cfdUpdateTransientTermCoefficients(theEquationName);
end
%
% update Property Fields
%
thePropertyNames = cfdGetPropertyNames;
theNumberOfProperties = length(thePropertyNames);
for iProperty = 1:theNumberOfProperties
    thePropertyName = thePropertyNames{iProperty};
    cfdUpdateTransientProperty(thePropertyName);
end

