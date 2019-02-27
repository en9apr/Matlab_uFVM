
function cfdUpdateFields
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%
%--------------------------------------------------------------------------

% update Property Fields
%------------------------------
% update BC for scalar Fields
%------------------------------
%
cfdUpdateProperties;
%
% update BC for scalar Fields
%
cfdUpdateEquations;
%
% Normalize volume fraction fields
%
applicationClass = cfdGetApplicationClass;
if strcmp(applicationClass, 'mutliphase')
    cfdNormalizeVFFields;
end