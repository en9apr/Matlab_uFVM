function cfdUpdateEquations
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   
%--------------------------------------------------------------------------
    
theEquationNames = cfdGetEquationNames;
theNumberOfEquations = length(theEquationNames);

for iEquation=1:theNumberOfEquations
    theEquationName = theEquationNames{iEquation};    

    cfdUpdateBoundaryConditions(theEquationName);
    
    cfdUpdateGradient(theEquationName);
    
    cfdUpdateTermCoefficientField(theEquationName)
    
    cfdUpdateScale(theEquationName);
end
