function theFluidNames = cfdGetAllFluidNames
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theNumberOfFluids = cfdGetNumberOfFluids;
theFluids = cfdGetFluids;

theFluidNames = {theFluids{1}.name};
for iFluid=2:theNumberOfFluids
    
theFluidNames = {theFluidNames theFluids{iFluid}.name};
    
end