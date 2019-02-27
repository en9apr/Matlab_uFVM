function theFluidName= cfdGetFluidNameUsingIndex(theFluidIndex)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theFluids = cfdGetFluids;

theFluid = theFluids{theFluidIndex};


theFluidName = theFluid.name;
