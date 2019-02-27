function cfdAssembleTransientTerm(theEquationName,theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function directs the solver towards the rigth time scheme
%--------------------------------------------------------------------------

% Get time scheme
theScheme = theTerm.scheme;

if strcmp(theScheme,'steadyState')
    return;
elseif strcmp(theScheme,'Euler')
    cfdAssembleTransientTermEuler(theEquationName,theTerm,iComponent);
elseif strcmp(theScheme,'Backward')
%     cfdAssembleTransientTermBackward(theEquationName,theTerm,iComponent);   
elseif strcmp(theScheme,'CrankNicolson')
%     cfdAssembleTransientTermCrankNicolson(theEquationName,theTerm,iComponent);
end