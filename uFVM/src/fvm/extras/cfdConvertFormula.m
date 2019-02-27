function theConvertedFormula = cfdConvertFormula(theFormula)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theConvertedFormula = theFormula;

theNumberOfFluids = cfdGetNumberOfFluids;
for iFluid=1:theNumberOfFluids
    theFluid = cfdGetFluidUsingIndex(iFluid);
    theFluidName = theFluid.name;
    theFluidTag = theFluid.tag;
    theConvertedFormula = regexprep(theConvertedFormula,[':' theFluidName],theFluidTag);
end


