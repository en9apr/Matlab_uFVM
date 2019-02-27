function  [initialResidual, finalResidual] = cfdSolveEquation(theEquationName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   
%--------------------------------------------------------------------------

theEquation = cfdGetModel(theEquationName);
isMultigridActive = theEquation.multigrid.isActive;
%
if isMultigridActive
    [initialResidual, finalResidual] = cfdApplyAMG(theEquation);
else
    [initialResidual, finalResidual] = cfdSolver(theEquation);
end
