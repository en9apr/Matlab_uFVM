function [rmsResidual,maxResidual,initialResidual,finalResidual] = cfdAssembleAndCorrectEquation(theEquationName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles and corrects the equation
%--------------------------------------------------------------------------

theEquation = cfdGetModel(theEquationName);
theEquationType = theEquation.type;

if strcmp(theEquationType,'Scalar')
    theNumberOfComponents = 1;
elseif strcmp(theEquationType,'Vector')
    theNumberOfComponents = 3;
end

% Get number of correctors
global Domain;
if strcmp(theEquationName, 'p')
    algorithm = cfdGetAlgorithm;
    nCorrectors = Domain.foam.fvSolution.(algorithm).nCorrectors;
else
    nCorrectors = 1;
end

% Initialize residuals of the equation
rmsResidual(1:theNumberOfComponents) = 0;
maxResidual(1:theNumberOfComponents) = 0;
initialResidual(1:theNumberOfComponents) = 0;
finalResidual(1:theNumberOfComponents) = 0;

for iComponent=1:theNumberOfComponents
    for iCorrector=1:nCorrectors
        % Assemble Equation
        [rmsResidual(iComponent), maxResidual(iComponent)] = cfdAssembleEquation(theEquationName,iComponent);
        
        % Solve Equation
        [initialResidual(iComponent),finalResidual(iComponent)] = cfdSolveEquation(theEquationName);
        
        % Correct Equation
        cfdCorrectEquation(theEquationName,iComponent);
    end
    
    % Store rms residual in the equation model
    cfdStoreResiduals(theEquationName, iComponent, rmsResidual(iComponent));    
end
