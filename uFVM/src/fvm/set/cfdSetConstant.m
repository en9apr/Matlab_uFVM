function cfdSetConstant(theConstant)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function stores the constant in the data base
%--------------------------------------------------------------------------

global Domain;

theIndex = 0;

theNumberOfConstants = length(Domain.constants);
for iConstant=1:theNumberOfConstants
    if strcmp(theConstant.name, Domain.constants{iConstant}.name)
        theIndex = iConstant;
    end
end

if theIndex==0
    theNumberOfConstants = length(Domain.constants);
    Domain.constants{theNumberOfConstants+1} = theConstant;
else
    disp('Overriding constant definition');
    Domain.constants{theIndex} = theConstant;
end
