function res = fvmAssembleAndCorrectEquationExplicitly(theScalarUserName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

% Assemble and Correct
%------------------------------
theScalarField = cfdGetModel(theScalarUserName);
theNumberOfTerms = length(theScalarField.terms);
if(theNumberOfTerms > 0)
    %Problem after the following function where Mdot is set to zero
    %*************************************************************    
    % Assemble Equation
    %------------------------------    
    res = fvmAssembleEquation(theScalarUserName);    
    %
    % Solve Equation
    %------------------------------
     % Solve Equation
    %------------------------------
    fvmSolveEquationExplicitly(theScalarUserName);
    %
    % Correct Equation
    %------------------------------
    fvmCorrectEquation(theScalarUserName);
    %
    % PPost Solve Updates
    %------------------------------
    fvmPostSolveEquation(theScalarUserName);
    %
    % Compute Residual
    %
    theResidual = fvmComputeResidual(theScalarUserName);
end
