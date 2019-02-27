function tolerance = cfdGetResidualTolerance(theEquationName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function retreives residual tolerance for each equation from
%   domain
%--------------------------------------------------------------------------

global Domain;
algorithm = cfdGetAlgorithm;
residualControl = Domain.foam.fvSolution.(algorithm).residualControl;
for iEquation=1:size(residualControl, 1)
   if strcmp(residualControl{iEquation, 1}, theEquationName)
       tolerance = str2double(residualControl{iEquation, 2});
       break;
   end
end