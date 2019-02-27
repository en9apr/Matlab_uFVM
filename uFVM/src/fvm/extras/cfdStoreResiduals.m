function cfdStoreResiduals(theEquationName, iComponent, rmsResidual)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function stores the residuals in the equation model
%--------------------------------------------------------------------------

% Get the current iteration from the length of the residuals array stored
% in the equation model
theEquation = cfdGetModel(theEquationName);
theSize = size(theEquation.residuals);
iteration = theSize(1);
if(iComponent==1)
    iteration = iteration + 1;
end

% Store the RMS residual in the equation model
theEquation.residuals(iteration,iComponent) = rmsResidual;
cfdSetModel(theEquation);