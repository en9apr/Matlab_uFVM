function theFluid = cfdGetFluidUsingName(theFluidName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theFluids = cfdGetFluids;
theNumberOfFluids = cfdGetNumberOfFluids;
theFluid = '';

for iFluid=1:theNumberOfFluids
    if strcmp(theFluidName, theFluids{iFluid}.name)
       theFluid = theFluids{iFluid};
    end
end