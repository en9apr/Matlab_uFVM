function theFluid = cfdGetFluid(theFluidName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theFluid={};
theFluids = cfdGetFluids;
theNumberOfFluids = cfdGetNumberOfFluids;
for iFluid=1:theNumberOfFluids
    if(strcmp(theFluidName,theFluids{iFluid}.name) || strcmp(theFluidName,theFluids{iFluid}.tag))
        theFluid = theFluids{iFluid};
        return
    end
end
