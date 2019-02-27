function theTag = cfdGetFluidTagUsingIndex(theFluidIndex)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================


theFluids = cfdGetFluids;

theFluid = theFluids{theFluidIndex};


theTag = theFluid.tag;
