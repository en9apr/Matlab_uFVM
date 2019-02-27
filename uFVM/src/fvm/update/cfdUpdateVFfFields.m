function cfdUpdateVFfFields
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   
%--------------------------------------------------------------------------

%
% Update VF_f field
%
theNumberOfFluids = cfdGetNumberOfFluids;
for iFluid=1:theNumberOfFluids 
    cfdUpdateVFfWithFluidIndex(iFluid);
end

