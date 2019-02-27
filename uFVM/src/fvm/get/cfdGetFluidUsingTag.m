function theFluid = cfdGetFluidUsingTag(theFluidTag)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theFluids = cfdGetFluids;

theNumberOfFluids = cfdGetNumberOfFluids;

for iFluid=1:theNumberOfFluids

    if(strcmp(theFluidTag,theFluids{iFluid}.tag))
    
theFluid = theFluids{iFluid};
    end
end