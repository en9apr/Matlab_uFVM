function theEquation = cfdAddTerm(theEquationName,theTermName,varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function sets up the terms of the equation and stores them
%--------------------------------------------------------------------------

theNumberOfSubTerms = length(varargin);

theEquationName = cfdConvertName(theEquationName);
theEquation = cfdGetModel(theEquationName);
theTerm = cfdSetupTerm(theTermName);

theTerm.fvmOperator = '';
theTerm.coefficientName = '';
theTerm.variableName = theEquationName;
theTerm.scheme = 'default';
if strcmp(theTermName,'Source')
    theTerm.formula = '';
end

for iSubTerm = 1:2:theNumberOfSubTerms
    theTerm =  setfield(theTerm,varargin{iSubTerm},varargin{iSubTerm+1});
end

if strcmp(theTermName,'Transient')
    theEquation.rhoName = theTerm.coefficientName;
end
if strcmp(theTermName,'Convection')
    theEquation.psiName = theTerm.coefficientName;
end
if strcmp(theTermName,'Diffusion') || strcmp(theTermName,'Stress')
    theEquation.gammaName = theTerm.coefficientName;
end

theNumberOfTerms = length(theEquation.terms);
theEquation.terms{theNumberOfTerms+1} =  theTerm ;

cfdSetModel(theEquation);