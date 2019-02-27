function theValue = cfdGetConstant(theConstantName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function retrieves the constant value from the data base
%--------------------------------------------------------------------------
global Domain;

theValue = 0;

theNumberOfConstants = length(Domain.constants);
for iConstant=1:theNumberOfConstants
    if strcmp(theConstantName,Domain.constants{iConstant}.name)
        theValue = Domain.constants{iConstant}.value;
        return;
    end
end

