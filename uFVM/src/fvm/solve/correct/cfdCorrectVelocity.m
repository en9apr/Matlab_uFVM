function cfdCorrectVelocity
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theNumberOfFluids = cfdGetNumberOfFluids;
for iFluid=1:theNumberOfFluids
    
    
    theFluid = cfdGetFluidUsingIndex(iFluid);
    theFluidTag = theFluid.tag;
  
    theMdotField_Before = cfdGetField(['Mdot' theFluidTag]);
    
    theMdotField = cfdCorrectMdotField(theFluidTag);
    cfdSetMeshField(theMdotField);
    
    save 'multi'  theMdotField_Before theMdotField
end
end
  