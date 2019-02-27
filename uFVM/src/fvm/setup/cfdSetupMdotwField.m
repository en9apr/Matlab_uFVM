function cfdSetupMdotwField

theNumberOfFluids = cfdGetNumberOfFluids;

for iFluid=1:theNumberOfFluids
    theFluidTag = cfdGetFluidTagUsingIndex(iFluid);
    theMdotName = ['Mdot_W' theFluidTag];
    cfdSetupMeshField(theMdotName,'Faces');
    %
    %VF_f = cfdSetupMeshField(['VF_f' theFluidTag],'Faces');
end
