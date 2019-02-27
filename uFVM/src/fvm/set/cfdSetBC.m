function theEquation = cfdSetBC(theEquationName,theBCPatch,varargin)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function sets up boundary conditions of the given equation
%--------------------------------------------------------------------------
theConvertedEquationName = cfdConvertName(theEquationName);
theEquation = cfdGetModel(theConvertedEquationName);

theBC.index = theBCPatch;
% default bc
if strcmp(theEquation.type,'Vector')
    theBC.type = 'noSlip';
    theBC.value = '[0;0;0]';
else
    theBC.type = 'specifiedValue';
    theBC.value = '0';
end

% setup the bc
theNumberOfSubTerms = length(varargin);
for iSubTerm = 1:2:theNumberOfSubTerms
    theBC = setfield(theBC,varargin{iSubTerm},varargin{iSubTerm+1});
end

% update equation and store back into database
theEquation.bcs{theBCPatch} = theBC;

cfdSetModel(theEquation);
