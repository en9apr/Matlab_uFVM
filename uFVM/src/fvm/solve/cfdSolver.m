function [initialResidual, finalResidual] = cfdSolver(theEquation)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   
%--------------------------------------------------------------------------

smootherType = theEquation.theSolutionMethods.smootherType;
[initialResidual, finalResidual] = cfdSolveAlgebraicSystem(1,smootherType,10);




