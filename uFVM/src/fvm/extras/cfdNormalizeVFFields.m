function cfdNormalizeVFFields
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   Normalize VF Fields at CV centers
%--------------------------------------------------------------------------
 
theNumberOfFluids = cfdGetNumberOfFluids;
SumVF = 0;
for iFluid=1:theNumberOfFluids
    VF = cfdGetMeshField(['alpha', num2str(iFluid)]);
    SumVF = SumVF + VF.phi;
end
for iFluid=1:theNumberOfFluids
    VF = cfdGetMeshField(['alpha', num2str(iFluid)]);
    VF.phi = VF.phi./SumVF;
    cfdSetMeshField(VF);
end


