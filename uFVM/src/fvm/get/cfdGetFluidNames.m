function theFluidNames = cfdGetFluidNames
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theNumberOfFluids = cfdGetNumberOfFluids;
theFluids = cfdGetFluids;

theFluidNames = {};
for iFluid=1:theNumberOfFluids    
    theFluidNames{iFluid} = theFluids{iFluid}.name;    
end